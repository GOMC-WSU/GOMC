/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef TARGETEDSWAP_H
#define TARGETEDSWAP_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "ConfigSetup.h"
#include <cmath>

//struct config_setup::TargetSwapParam;

struct TSwapParam {
  // defines the center of subVolume
  XYZ subVolumeCenter; 
  // defines the dimension of subVolume
  XYZ subVolumeDim; 
  // defines the targeted molecule kind
  std::vector<uint> selectedResKind; 
  // defines the subVolume index for error checking
  int subVolumeIdx; 
  // volume of subVolume (A^3)
  double subVolume;
};

class TargetedSwap : public MoveBase
{
public:

  TargetedSwap(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
    MoveBase(sys, statV)
    {
      for(int b = 0; b < BOX_TOTAL; b++) {
        hasSubVolume[b] = false;
        pickedSubV[b] = 0;
      }

      if(statV.targetedSwapVal.enable) {
        for(int s = 0; s < statV.targetedSwapVal.targetedSwap.size(); s++) {
          config_setup::TargetSwapParam tsp = statV.targetedSwapVal.targetedSwap[s];
          TSwapParam tempVar;
          // copy data from TargetSwapParam struct to TSwapParam
          tempVar.subVolumeCenter = tsp.subVolumeCenter;
          tempVar.subVolumeDim = tsp.subVolumeDim;
          tempVar.subVolumeIdx = tsp.subVolumeIdx;
          tempVar.subVolume = tsp.subVolumeDim.x * tsp.subVolumeDim.y *
                              tsp.subVolumeDim.z;
          // Use all residue kind
          if (tsp.selectedResKind[0] == "ALL") {
            for(uint k = 0; k < molLookRef.GetNumCanSwapKind(); k++) {
              uint kindIdx = molLookRef.GetCanSwapKind(k);
              tempVar.selectedResKind.push_back(kindIdx);
            }
          } else {
            // Sort user defined residue kind
            std::sort(tsp.selectedResKind.begin(), tsp.selectedResKind.end());
            // Loop through all user defined residue name
            for (uint i = 0; i < tsp.selectedResKind.size(); i++) {
              std::string kindName = tsp.selectedResKind[i];
              bool found = false;
              // Check if we have "kindName" in system
              for(uint k = 0; k < molLookRef.GetNumCanSwapKind(); k++) {
                uint kindIdx = molLookRef.GetCanSwapKind(k);
                if(molRef.kinds[kindIdx].name == kindName) {
                  found = true;
                  tempVar.selectedResKind.push_back(kindIdx);
                  break; // break from first loop
                }
              }
              // If there was no such a residue name, through error
              if(!found) {
                printf("Error: In Targeted Swap move, residue name %s cannot be swapped or was not found in PDB/PSF file!\n",
                        kindName.c_str());
                exit(EXIT_FAILURE);
              }
            }
          }
          targetSwapParam[tsp.selectedBox].push_back(tempVar);
          hasSubVolume[tsp.selectedBox] = true;
        }
      }
    }

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);
  virtual void PrintAcceptKind();

private:

  double GetCoeff() const;
  uint GetBoxPairAndMol(const double subDraw, const double movPerc);
  uint PickMolInSubVolume(const uint &box);
  void CheckSubVolumePBC(const uint &box);
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;

  std::vector<TSwapParam> targetSwapParam[BOX_TOTAL];
  std::vector<uint> molIdxInSubVolume[BOX_TOTAL];
  XYZ subVolCenter[BOX_TOTAL], subVolDim[BOX_TOTAL];
  double subVolume[BOX_TOTAL];
  bool hasSubVolume[BOX_TOTAL];
  int pickedSubV[BOX_TOTAL];
};

void TargetedSwap::PrintAcceptKind()
{
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Targeted-Swap ", molRef.kinds[k].name.c_str());
    for(uint b = 0; b < BOX_TOTAL; b++) {
      if(moveSetRef.GetTrial(b, mv::TARGETED_SWAP, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::TARGETED_SWAP, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint TargetedSwap::Prep(const double subDraw, const double movPerc)
{
  overlap = false;
  // Pick the source, dest box. Pick the subVolume in each box (if exist)
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStart);

    // We need to copy the oldMol coordinate to newMol in order to be able to
    // insert rigidly in dest box. This means we need to unwrap it using source box,
    // and then wrap it using dest box
    XYZArray molA(pLen);
    coordCurrRef.CopyRange(molA, pStart, 0, pLen);
    boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndex));
    boxDimRef.WrapPBC(molA, destBox);
    //set coordinate of mole to newMol, later it will shift to center
    newMol.SetCoords(molA, 0);

    // the trial configuration has cavity but COM is not fix and no rotation
    // around backbone
    if(hasSubVolume[sourceBox]) {
      // Since the subVolume is orthogonal, no need to calculate and set the cellBasis
      // matrix. The default value for cellBasis matrix is I 
      oldMol.SetSeed(subVolCenter[sourceBox], subVolDim[sourceBox], true, false, false);
    }
    if(hasSubVolume[destBox]) {
      // Since the subVolume is orthogonal, no need to calculate and set the cellBasis
      // matrix. The default value for cellBasis matrix is I 
      newMol.SetSeed(subVolCenter[destBox], subVolDim[destBox], true, false, false);
    }
  }
  return state;
}

inline uint TargetedSwap::GetBoxPairAndMol(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;

  // 1. Pick a pair of box to determin source and destination box.
  prng.PickBoxPair(sourceBox, destBox, subDraw, movPerc);
  molIdxInSubVolume[sourceBox].clear();
  molIdxInSubVolume[destBox].clear();
  // If we have a subvolume for source box, pick one of the subvolume randomely
  if(hasSubVolume[sourceBox]) {
    int SubVIdx_source = prng.randIntExc(targetSwapParam[sourceBox].size());
    // Set the picked subVolume parameter for source box
    pickedSubV[sourceBox] = SubVIdx_source;
    subVolCenter[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolumeCenter;
    subVolDim[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolumeDim;
    subVolume[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolume;
    // Wrap the center of subVolume and make sure the it's dimension is less than
    // simulation box
    CheckSubVolumePBC(sourceBox);
  }

  // If we have a subvolume for dest box, pick one of the subvolume randomely
  if(hasSubVolume[destBox]) {
    int SubVIdx_dest = prng.randIntExc(targetSwapParam[destBox].size());
    // Set the picked subVolume parameter for dest box
    pickedSubV[destBox] = SubVIdx_dest;
    subVolCenter[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumeCenter;
    subVolDim[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumeDim;
    subVolume[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolume;
    // Wrap the center of subVolume and make sure the it's dimension is less than
    // simulation box
    CheckSubVolumePBC(destBox);
  }

  // 2. Pick a molecule kind
  // For the simulations, where all subVolumes are in one box we can do the following.
  // If there are subVolumes in each box, then we have to make sure that targeted 
  // residue kind in all subvolumes are exactly identical. Otherwise the results would
  //  would be biased.
  if(hasSubVolume[sourceBox]) {
    // Pick a molecule kind that specified for sourceBox subvolume
    uint i = prng.randIntExc(targetSwapParam[sourceBox][pickedSubV[sourceBox]].selectedResKind.size());
    kindIndex = targetSwapParam[sourceBox][pickedSubV[sourceBox]].selectedResKind[i];
  } else if (hasSubVolume[destBox]) {
    // Pick a molecule kind that specified for destBox subvolume
    uint i = prng.randIntExc(targetSwapParam[destBox][pickedSubV[destBox]].selectedResKind.size());
    kindIndex = targetSwapParam[destBox][pickedSubV[destBox]].selectedResKind[i];
  } else {
    printf("Error: In Targeted Swap move, no subVolume was defined for any Box!\n");
    exit(EXIT_FAILURE);
  }

  // 3. Pick a molecule Index for the picked molecule kind
  if(hasSubVolume[sourceBox]) {
    // Search through all molecule in the subVolume and randomely pick 
    // one of the picked kind
    state = PickMolInSubVolume(sourceBox);
  } else {
    // Randomely pick one molecule of the picked kind from whole source box
    state = prng.PickMolIndex(molIndex, kindIndex, sourceBox);
    // Just call the function to get number of molecule in cavity in destBox
    calcEnRef.FindMolInCavity(molIdxInSubVolume[destBox],
                              subVolCenter[destBox], subVolDim[destBox],
                              destBox, kindIndex);
  }

#if ENSEMBLE == GCMC
  if(state == mv::fail_state::NO_MOL_OF_KIND_IN_BOX && sourceBox == mv::BOX1) {
    std::cout << "Error: There are no molecules of kind " <<
              molRef.kinds[kindIndex].name << " left in reservoir.\n";
    exit(EXIT_FAILURE);
  }
#endif

  // If we found the molecule, set the start and length of the molecule
  if (state == mv::fail_state::NO_FAIL) {
    pStart = pLen = 0;
    molRef.GetRangeStartLength(pStart, pLen, molIndex);
  }
  return state;
}

inline void TargetedSwap::CheckSubVolumePBC(const uint &box)
{
  XYZ boxAxis = boxDimRef.GetAxis(box);
  bool shrunk = false;
  // Check the subVolume dimension
  if(subVolDim[box].x > boxAxis.x) {
    printf("Error: In Targeted Swap move, box dimension in X axis shrunk bellow %g A for Subvolume %d!\n",
           subVolDim[box].x, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }
  if(subVolDim[box].y > boxAxis.y) {
    printf("Error: In Targeted Swap move, box dimension in Y axis shrunk bellow %g A for Subvolume %d!\n",
           subVolDim[box].y, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }
  if(subVolDim[box].z > boxAxis.z) {
    printf("Error: In Targeted Swap move, box dimension in Z axis shrunk bellow %g A for Subvolume %d!\n",
           subVolDim[box].z, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }

  // Wrap the center of subVolume and check for PBC
  subVolCenter[box] = boxDimRef.WrapPBC(subVolCenter[box], box);
  // Check the center or subvolume
  if(subVolCenter[box].x >= boxAxis.x) {
    printf("Error: In Targeted Swap move, box dimension in X axis shrunk bellow center %g A for Subvolume %d!\n",
           subVolCenter[box].x, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }
  if(subVolCenter[box].y >= boxAxis.y) {
    printf("Error: In Targeted Swap move, box dimension in Y axis shrunk bellow center %g A for Subvolume %d!\n",
           subVolCenter[box].y, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }
  if(subVolCenter[box].z >= boxAxis.z) {
    printf("Error: In Targeted Swap move, box dimension in Z axis shrunk bellow center %g A for Subvolume %d!\n",
           subVolCenter[box].z, targetSwapParam[box][pickedSubV[box]].subVolumeIdx);
    shrunk = true;
  }

  if (shrunk) {
    exit(EXIT_FAILURE);
  }
}

inline uint TargetedSwap::PickMolInSubVolume(const uint &box)
{
  uint rejectState = mv::fail_state::NO_FAIL;
  // Find all the molecule of kindIndex in subvolume
  bool foundMol = calcEnRef.FindMolInCavity(molIdxInSubVolume[box],
                  subVolCenter[box], subVolDim[box], box, kindIndex);
  if(foundMol) {
    // The return vector, stores unique molecule index
    uint i = prng.randIntExc(molIdxInSubVolume[box].size());
    molIndex = molIdxInSubVolume[box][i];
  } else {
    //reject the move if no molecule was find
    rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
  }
  
  return rejectState;
}

inline uint TargetedSwap::Transform()
{
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  molRef.kinds[kindIndex].BuildIDOld(oldMol, molIndex);
  //Add bonded energy because we don't consider it in DCRotate.cpp
  oldMol.AddEnergy(calcEnRef.MoleculeIntra(oldMol));

  molRef.kinds[kindIndex].BuildIDNew(newMol, molIndex);
  //Add bonded energy because we don't consider it in DCRotate.cpp
  newMol.AddEnergy(calcEnRef.MoleculeIntra(newMol));

  overlap = newMol.HasOverlap();
  return mv::fail_state::NO_FAIL;
}

inline void TargetedSwap::CalcEn()
{
  W_tc = 1.0;
  W_recip = 1.0;
  correct_old = 0.0;
  correct_new = 0.0;
  self_old = 0.0;
  self_new = 0.0;

  if (ffRef.useLRC) {
    tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
    tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
    W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
  }

  if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
    correct_new = calcEwald->SwapCorrection(newMol);
    correct_old = calcEwald->SwapCorrection(oldMol);
    self_new = calcEwald->SwapSelf(newMol);
    self_old = calcEwald->SwapSelf(oldMol);
    //SwapDestRecip must be called first to backup the cosMol and sinMol
    recipGain.energy =
      calcEwald->SwapDestRecip(newMol, destBox, molIndex);
    recipLose.energy =
      calcEwald->SwapSourceRecip(oldMol, sourceBox, molIndex);
    //need to contribute the self and correction energy
    W_recip = exp(-1.0 * ffRef.beta * (recipGain.energy + recipLose.energy +
                                       correct_new - correct_old +
                                       self_new - self_old));
  }

}

inline double TargetedSwap::GetCoeff() const
{
#if ENSEMBLE == GEMC
  printf("Error: GOMC does not support Targeted swap for GEMC!\n");
  exit(EXIT_FAILURE);
  return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) /
         (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
         boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) { //Delete case
    double molNumber =  double(molIdxInSubVolume[sourceBox].size());
    double invVolume = 1.0 / subVolume[sourceBox];
    /*
    if (molIdxInSubVolume[sourceBox].size() != molLookRef.NumKindInBox(kindIndex, sourceBox)) {
      printf("Error: %d in subVolume vs %d in box!\n", molIdxInSubVolume[sourceBox].size(),
      molLookRef.NumKindInBox(kindIndex, sourceBox));
    } */
    if(ffRef.isFugacity) {
      return molNumber * invVolume /
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return molNumber * invVolume *
             exp(-BETA * molRef.kinds[kindIndex].chemPot);
    }
  } else { //Insertion case
    double molNumber =  double(molIdxInSubVolume[destBox].size());
    double volume = subVolume[destBox];
    /*
    if (molIdxInSubVolume[destBox].size() != molLookRef.NumKindInBox(kindIndex, destBox)) {
      printf("Error: %d in subVolume vs %d in box!\n", molIdxInSubVolume[destBox].size(),
      molLookRef.NumKindInBox(kindIndex, destBox));
    } */
    if(ffRef.isFugacity) {
      return volume / (molNumber + 1) *
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return volume / (molNumber + 1) *
             exp(BETA * molRef.kinds[kindIndex].chemPot);
    }
  }
#endif
}

inline void TargetedSwap::Accept(const uint rejectState, const uint step)
{
  bool result;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    double molTransCoeff = GetCoeff();
    double Wo = oldMol.GetWeight();
    double Wn = newMol.GetWeight();
    double Wrat = Wn / Wo * W_tc * W_recip;

    //safety to make sure move will be rejected in overlap case
    if(newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
      result = prng() < molTransCoeff * Wrat;
    } else
      result = false;

    if(result) {
      //Add tail corrections
      sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
      sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
      //Add rest of energy.
      sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
      sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
      //Add Reciprocal energy
      sysPotRef.boxEnergy[sourceBox].recip += recipLose.energy;
      sysPotRef.boxEnergy[destBox].recip += recipGain.energy;
      //Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
      sysPotRef.boxEnergy[destBox].correction += correct_new;
      //Add self energy
      sysPotRef.boxEnergy[sourceBox].self -= self_old;
      sysPotRef.boxEnergy[destBox].self += self_new;

      //Set coordinates, new COM; shift index to new box's list
      newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
      comCurrRef.SetNew(molIndex, destBox);
      molLookRef.ShiftMolBox(molIndex, sourceBox, destBox,
                             kindIndex);
      cellList.AddMol(molIndex, destBox, coordCurrRef);

      //Zero out box energies to prevent small number
      //errors in double.
      if (molLookRef.NumInBox(sourceBox) == 0) {
        sysPotRef.boxEnergy[sourceBox].Zero();
        sysPotRef.boxVirial[sourceBox].Zero();
      } else if (molLookRef.NumInBox(sourceBox) == 1) {
        sysPotRef.boxEnergy[sourceBox].inter = 0;
        sysPotRef.boxVirial[sourceBox].inter = 0;
        sysPotRef.boxEnergy[sourceBox].real = 0;
        sysPotRef.boxVirial[sourceBox].real = 0;
      }

      calcEwald->UpdateRecip(sourceBox);
      calcEwald->UpdateRecip(destBox);

      //Retotal
      sysPotRef.Total();
    } else {
      cellList.AddMol(molIndex, sourceBox, coordCurrRef);
      //when weight is 0, MolDestSwap() will not be executed, thus cos/sin
      //molRef will not be changed. Also since no memcpy, doing restore
      //results in memory overwrite
      if (newMol.GetWeight() != 0.0 && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }

  } else //we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::TARGETED_SWAP, result, destBox, kindIndex);
}

#endif

#endif
