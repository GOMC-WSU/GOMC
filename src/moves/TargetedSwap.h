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
#include "Geometry.h" // for bond info
#include "FloydWarshallCycle.h"
#include <cmath>
#include <queue>

//struct config_setup::TargetSwapParam;
struct BondList;

struct TSwapParam {
  // defines the center of subVolume
  XYZ subVolumeCenter; 
  // defines the dimension of subVolume
  XYZ subVolumeDim; 
  // defines the targeted molecule kind
  std::vector<uint> selectedResKind; 
  // defines the list if atom to calculcate the center of
  // subVolume
  std::vector<int> atomList; 
  // defines the subVolume index for error checking
  int subVolumeIdx; 
  // volume of subVolume (A^3)
  double subVolume;
  // swap type rigid/flexible
  bool rigidSwap;
  // use atomList to find subvolume center
  bool calcSubVolCenter;
};

class TargetedSwap : public MoveBase
{
public:

  TargetedSwap(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
    MoveBase(sys, statV)
    {
      rigidSwap = true;
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
          tempVar.atomList = tsp.atomList;
          tempVar.subVolumeDim = tsp.subVolumeDim;
          tempVar.subVolumeIdx = tsp.subVolumeIdx;
          tempVar.subVolume = tsp.subVolumeDim.x * tsp.subVolumeDim.y *
                              tsp.subVolumeDim.z;
          tempVar.rigidSwap = tsp.rigid_swap;
          // if the size of atom list is not zero, then we must calculate the center
          // of subVolume
          tempVar.calcSubVolCenter = tempVar.atomList.size();
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

        //Initialize the trial and acceptance
        for(uint b = 0; b < BOX_TOTAL; ++b) {
          if(hasSubVolume[b]) {
            trial[b].resize(targetSwapParam[b].size());
            accepted[b].resize(targetSwapParam[b].size());
            for(uint idx = 0; idx < targetSwapParam[b].size(); ++idx) {
              trial[b][idx].resize(molRef.GetKindsCount(), 0);
              accepted[b][idx].resize(molRef.GetKindsCount(), 0);
            }
          }
        }

        // Get the atom index for growing seed for each atomKind
        for(uint k = 0; k < molLookRef.GetNumKind(); k++) {
          growingAtomIndex.push_back(GetGrowingAtomIndex(k));
        }

        PrintTargetedSwapInfo();
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
  uint PickMolInSubVolume();
  // Calculate the subvolume center, using defined atomList
  XYZ GetSubVolumeCenter(const uint &box, const uint &SubVIdx);
  // Check periodic boundary for subVolume
  void CheckSubVolumePBC(const uint &box);
  void  PrintTargetedSwapInfo();
  // Track the acceptance for each subVolume
  void AcceptSubVolumeIdx(const uint &state, const uint &kind, const uint &box);
  // Returns the center node in molecule's graph
  uint GetGrowingAtomIndex(const uint k);
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;

  std::vector<uint> growingAtomIndex;
  std::vector<TSwapParam> targetSwapParam[BOX_TOTAL];
  std::vector<uint> molIdxInSubVolume[BOX_TOTAL];
  //For move acceptance of each molecule subVolume, and kind
  std::vector< std::vector<uint> > trial[BOX_TOTAL], accepted[BOX_TOTAL];
  XYZ subVolCenter[BOX_TOTAL], subVolDim[BOX_TOTAL];
  double subVolume[BOX_TOTAL];
  bool hasSubVolume[BOX_TOTAL];
  int pickedSubV[BOX_TOTAL];
  bool rigidSwap;
};

// Need to fix it for GEMC
void TargetedSwap::PrintAcceptKind()
{
  for(uint b = 0; b < BOX_TOTAL; ++b) {
    if(hasSubVolume[b]) {
      for(uint idx = 0; idx < targetSwapParam[b].size(); ++idx) {
        for(uint k = 0; k < targetSwapParam[b][idx].selectedResKind.size(); ++k) {
          int t = trial[b][idx][k];
          int a = accepted[b][idx][k];
          double percentAccept = (t ? 100 * (double)a/t: 0.0);
          char msg[257];
          sprintf(msg, "%% Accepted Targeted-Swap (%d)", idx);
          printf("%-30s %-5s ", msg, molRef.kinds[k].name.c_str());
          printf("%10.5f %10.5f \n", percentAccept, percentAccept);
        }
      }
    }
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
    if(hasSubVolume[sourceBox]) { // deletion
      // Since the subVolume is orthogonal, no need to calculate and set the cellBasis
      // matrix. The default value for cellBasis matrix is I 
      oldMol.SetSeed(subVolCenter[sourceBox], subVolDim[sourceBox], true, false, false);
      oldMol.SetGrowingAtomIndex(growingAtomIndex[kindIndex]);
    }
    if(hasSubVolume[destBox]) { //insertion
      // Since the subVolume is orthogonal, no need to calculate and set the cellBasis
      // matrix. The default value for cellBasis matrix is I 
      newMol.SetSeed(subVolCenter[destBox], subVolDim[destBox], true, false, false);
      newMol.SetGrowingAtomIndex(growingAtomIndex[kindIndex]);
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
  rigidSwap = true;
  // If we have a subvolume for source box, pick one of the subvolume randomely
  if(hasSubVolume[sourceBox]) {
    int SubVIdx_source = prng.randIntExc(targetSwapParam[sourceBox].size());
    // Set the picked subVolume parameter for source box
    pickedSubV[sourceBox] = SubVIdx_source;
    subVolDim[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolumeDim;
    subVolume[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolume;
    rigidSwap = targetSwapParam[sourceBox][SubVIdx_source].rigidSwap;
    //determin if we should calculate the center of subvolume or not
    if (targetSwapParam[sourceBox][SubVIdx_source].calcSubVolCenter) {
      subVolCenter[sourceBox] = GetSubVolumeCenter(sourceBox, SubVIdx_source);
    } else {
      subVolCenter[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolumeCenter;
    }
    // Wrap the center of subVolume and make sure the it's dimension is less than
    // simulation box
    CheckSubVolumePBC(sourceBox);
  }

  // If we have a subvolume for dest box, pick one of the subvolume randomely
  if(hasSubVolume[destBox]) {
    int SubVIdx_dest = prng.randIntExc(targetSwapParam[destBox].size());
    // Set the picked subVolume parameter for dest box
    pickedSubV[destBox] = SubVIdx_dest;
    subVolDim[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumeDim;
    subVolume[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolume;
    rigidSwap = targetSwapParam[destBox][SubVIdx_dest].rigidSwap;
    //determin if we should calculate the center of subvolume or not
    if (targetSwapParam[destBox][SubVIdx_dest].calcSubVolCenter) {
      subVolCenter[destBox] = GetSubVolumeCenter(destBox, SubVIdx_dest);
    } else {
      subVolCenter[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumeCenter;
    }
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
  state = PickMolInSubVolume();
  

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

inline XYZ TargetedSwap::GetSubVolumeCenter(const uint &box, const uint &subVIdx)
{
  int i;
  TSwapParam tsp = targetSwapParam[box][subVIdx];
  int listSize = tsp.atomList.size();
  // use first atomIndex to unwrap others
  XYZ reference = coordCurrRef.Get(tsp.atomList[0]);
  XYZ center = reference;
  XYZ atomCoord;
  // start from next atom index
  for(i = 1; i < listSize; ++i) {
    if(molLookRef.IsAtomInBox(tsp.atomList[i], box)) {
      atomCoord = coordCurrRef.Get(tsp.atomList[i]);
      center += boxDimRef.UnwrapPBC(atomCoord, box, reference);
    } else {
      printf("Error: In Targeted Swap move, atom index %d does not exist in box %d!\n",
            tsp.atomList[i], box);
      printf("Error: Cannot calculate the center of subVolume %d!\n", subVIdx);
      exit(EXIT_FAILURE);
    }
  }
  center *= 1.0 / double(listSize);
  // wrap the center
  center = boxDimRef.WrapPBC(center, box);
  return center;
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

inline uint TargetedSwap::PickMolInSubVolume()
{

  uint rejectState = mv::fail_state::NO_FAIL;
  if(hasSubVolume[sourceBox]) {
    // Search through all molecule in the subVolume and randomely pick 
    // one of the picked kind
    bool foundMol = false;
    if(rigidSwap) {
      // count molecule in subVolume by molecule geometric center
      // Find all the molecule of kindIndex in subvolume
      foundMol = calcEnRef.FindMolInCavity(molIdxInSubVolume[sourceBox],
                    subVolCenter[sourceBox], subVolDim[sourceBox],
                    sourceBox, kindIndex);
    } else {
      // count molecule in subVolume by coordinate of seed atomIdx
      // Find all the molecule of kindIndex in subvolume
      foundMol = calcEnRef.FindMolInCavity(molIdxInSubVolume[sourceBox],
                    subVolCenter[sourceBox], subVolDim[sourceBox],
                    sourceBox, kindIndex, growingAtomIndex[kindIndex]);
    }

    if(foundMol) {
      // The return vector, stores unique molecule index
      uint i = prng.randIntExc(molIdxInSubVolume[sourceBox].size());
      molIndex = molIdxInSubVolume[sourceBox][i];
    } else {
      //reject the move if no molecule was find
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    }
  } else { // means we have subVolume in destBox
    // Randomely pick one molecule of the picked kind from whole source box
    rejectState = prng.PickMolIndex(molIndex, kindIndex, sourceBox);
    // Just call the function to get number of molecule in cavity in destBox
    if(rigidSwap) {
      // count molecule in subVolume by molecule geometric center
      // Find all the molecule of kindIndex in subvolume
      calcEnRef.FindMolInCavity(molIdxInSubVolume[destBox],
                subVolCenter[destBox], subVolDim[destBox],
                destBox, kindIndex);
    } else {
      // count molecule in subVolume by coordinate of seed atomIdx
      // Find all the molecule of kindIndex in subvolume
      calcEnRef.FindMolInCavity(molIdxInSubVolume[destBox],
                subVolCenter[destBox], subVolDim[destBox],
                destBox, kindIndex, growingAtomIndex[kindIndex]);
    }
  }

  return rejectState;
}

inline uint TargetedSwap::Transform()
{
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  if(rigidSwap) {
    molRef.kinds[kindIndex].BuildIDOld(oldMol, molIndex);
    //Add bonded energy because we don't consider it in DCRotateCOM.cpp
    oldMol.AddEnergy(calcEnRef.MoleculeIntra(oldMol));

    molRef.kinds[kindIndex].BuildIDNew(newMol, molIndex);
    //Add bonded energy because we don't consider it in DCRotateCOM.cpp
    newMol.AddEnergy(calcEnRef.MoleculeIntra(newMol));
  } else {
    // insert or delete using CD-CBMC starting from growingatom Index
    molRef.kinds[kindIndex].BuildGrowInCav(oldMol, newMol, molIndex);
  }
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
  AcceptSubVolumeIdx(result, kindIndex, destBox);
}


uint TargetedSwap::GetGrowingAtomIndex(const uint k)
{
  // if molecule kind is not swapable, dont waste time to get the
  // center node. This kind will not be selected to swap
  if (!molLookRef.IsKindSwapable(k)) {
    return 0;
  }

  // If molecule is mono or diatomic, retrune index 0
  if (molRef.kinds[k].NumAtoms() <= 2) { 
    return 0;
  }

  // need to find the center node.
  // https://codeforces.com/blog/entry/17974
  {
    MoleculeKind &kind = molRef.kinds[k];
    FloydWarshallCycle flw(kind.NumAtoms());
    //Setup the node's degree
    uint bondCount = kind.bondList.count;
    for (int i = 0; i < bondCount; ++i) {
      flw.AddEdge(kind.bondList.part1[i], kind.bondList.part2[i]);
    }

    // get the centric node
    int centric = flw.GetCentricNode();
    
    // If it was not successfull, it returns -1
    if(centric == -1) {
      printf("Error: In TargetedSwap move, no center atom was found for %s!\n",
              kind.name.c_str());
      exit(EXIT_FAILURE);
    }

    return centric;
  }
}

void TargetedSwap::PrintTargetedSwapInfo()
{
  int i, k, b;
  for(b = 0; b < BOX_TOTAL; ++b) {
    for (i = 0; i < targetSwapParam[b].size(); ++i) {
      TSwapParam tsp = targetSwapParam[b][i];
      printf("%-40s %d: \n",       "Info: Targeted Swap parameter for subVolume index",
              tsp.subVolumeIdx);
      printf("%-40s %d \n",       "      SubVolume Box:", b);
      if (tsp.calcSubVolCenter) {
        printf("%-40s Using %d defined atom indexe/s \n", "      Calculating subVolume center:",
                tsp.atomList.size());
        int max = *std::max_element(tsp.atomList.begin(), tsp.atomList.end());
        if(max >= coordCurrRef.Count()) {
          printf("Error: Atom index %d is beyond total number of atoms (%d) in the system!\n",
                 max, coordCurrRef.Count());
          printf("       Make sure to use 0 based atom index!\n");
          exit(EXIT_FAILURE);
        }
      } else {
        printf("%-40s %g %g %g \n", "      SubVolume center:",
                tsp.subVolumeCenter.x, tsp.subVolumeCenter.y, tsp.subVolumeCenter.z);
      }
      printf("%-40s %g %g %g \n", "      SubVolume dimension:",
              tsp.subVolumeDim.x, tsp.subVolumeDim.y, tsp.subVolumeDim.z);
      printf("%-40s %s \n",       "      SubVolume Swap type:", (tsp.rigidSwap ? "Rigid body" : "Flexible"));
      printf("%-40s ",            "      Targeted residue kind:");
      for (k = 0; k < tsp.selectedResKind.size(); ++k) {
        printf("%-5s ", molRef.kinds[k].name.c_str());
      }
      printf("\n");
      if(!tsp.rigidSwap) {
        for (k = 0; k < tsp.selectedResKind.size(); ++k) {
          printf("%-40s %s: (%d, %s) \n",       "      Starting atom (index, name) for",
                molRef.kinds[k].name.c_str(), growingAtomIndex[k],
                molRef.kinds[k].atomNames[growingAtomIndex[k]].c_str());
        }
      }
      printf("\n\n");
    }
  }
}


  void TargetedSwap::AcceptSubVolumeIdx(const uint &state, const uint &kind,
                                        const uint &box)
  {
    if (hasSubVolume[box]) {
      uint subVidx = pickedSubV[box];
      ++trial[box][subVidx][kind];
      if(state) {
        ++accepted[box][subVidx][kind];
      }
    }
  }

#endif

#endif
