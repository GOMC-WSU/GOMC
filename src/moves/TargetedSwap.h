/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef TARGETEDSWAP_H
#define TARGETEDSWAP_H

// enumerator for periodicity of each axis
enum PBC { PBC_X, PBC_Y, PBC_Z, PBC_XY, PBC_XZ, PBC_YZ, PBC_XYZ };

struct TSwapParam {
  // defines the center of subVolume
  XYZ subVolumeCenter;
  // defines the dimension of subVolume
  XYZ subVolumeDim;
  // defines the targeted molecule kind
  std::vector<uint> selectedResKind;
  // sub-Volume periodicity for each axis
  PBC subVolumePBC;
  // defines the list if atom to calculcate the center of subVolume
  std::vector<int> atomList;
  // defines the subVolume index for error checking
  int subVolumeIdx;
  // volume of subVolume (A^3)
  double subVolume;
  // swap type rigid/flexible
  bool rigidSwap;
  // use atomList to find subvolume center
  bool calcSubVolCenter;
  // sub-Volume periodicity for each axis
  std::vector<bool> xyzPBC;
#if ENSEMBLE == GCMC
  // to keep track of chemical potentials for each subVolume
  std::vector<double> chemPot;
#endif
};

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC

#include <queue>

#include "ConfigSetup.h"
#include "FloydWarshallCycle.h"
#include "GeomLib.h"
#include "Geometry.h" // for bond info
#include "MoveBase.h"
#include "TrialMol.h"

struct BondList;

class TargetedSwap : public MoveBase {
public:
  TargetedSwap(System &sys, StaticVals const &statV)
      : MoveBase(sys, statV), molLookRef(sys.molLookupRef),
        ffRef(statV.forcefield) {
    rigidSwap = true;
    for (int b = 0; b < BOX_TOTAL; b++) {
      hasSubVolume[b] = false;
      pbcMode[b] = PBC_XYZ;
      pickedSubV[b] = 0;
    }

    if (statV.targetedSwapVal.enable) {
#if ENSEMBLE == GEMC
      printf("Error: GOMC does not support Targeted-Swap for GEMC!\n");
      exit(EXIT_FAILURE);
#endif
      for (int s = 0; s < (int)statV.targetedSwapVal.targetedSwap.size(); s++) {
        config_setup::TargetSwapParam tsp =
            statV.targetedSwapVal.targetedSwap[s];
        TSwapParam tempVar;
        // copy data from TargetSwapParam struct to TSwapParam
        tempVar.subVolumeCenter = tsp.subVolumeCenter;
        tempVar.atomList = tsp.atomList;
        tempVar.subVolumeDim = tsp.subVolumeDim;
        tempVar.subVolumeIdx = tsp.subVolumeIdx;
        tempVar.subVolume =
            tsp.subVolumeDim.x * tsp.subVolumeDim.y * tsp.subVolumeDim.z;
        tempVar.rigidSwap = tsp.rigid_swap;

        // find and Shift the atom indicies for box > 0
        int ShiftAtomIndex = 0;
        for (int b = 0; b < tsp.selectedBox; b++) {
          for (int k = 0; k < (int)molRef.GetKindsCount(); k++) {
            MoleculeKind const &thisKind = molRef.kinds[k];
            ShiftAtomIndex +=
                (int)thisKind.NumAtoms() * sys.molLookupRef.NumKindInBox(k, b);
          }
        }
        for (int i = 0; i < tempVar.atomList.size(); i++) {
          tempVar.atomList[i] += ShiftAtomIndex;
        }

        // if the size of atom list is not zero, then we must calculate the
        // center of subVolume
        tempVar.calcSubVolCenter = tempVar.atomList.size();
        // Set the PBC info
        tempVar.xyzPBC = tsp.subVolumeBPC;
        if (tempVar.xyzPBC[0]) { // pbc in x axis
          tempVar.subVolumePBC = PBC_X;
          if (tempVar.xyzPBC[1]) { // pbc in y axis
            tempVar.subVolumePBC = PBC_XY;
            if (tempVar.xyzPBC[2]) { // pbc in z axis
              tempVar.subVolumePBC = PBC_XYZ;
            }
          } else if (tempVar.xyzPBC[2]) { // pbc in z axis
            tempVar.subVolumePBC = PBC_XZ;
          }
        } else if (tempVar.xyzPBC[1]) { // pbc in y axis
          tempVar.subVolumePBC = PBC_Y;
          if (tempVar.xyzPBC[2]) { // pbc in z axis
            tempVar.subVolumePBC = PBC_YZ;
          }
        } else if (tempVar.xyzPBC[2]) { // pbc in z axis
          tempVar.subVolumePBC = PBC_Z;
        }

        // Use all residue kind
        if (tsp.selectedResKind[0] == "ALL") {
          for (uint k = 0; k < molLookRef.GetNumCanSwapKind(); k++) {
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
            for (uint k = 0; k < molLookRef.GetNumCanSwapKind(); k++) {
              uint kindIdx = molLookRef.GetCanSwapKind(k);
              if (molRef.kinds[kindIdx].name == kindName) {
                found = true;
                tempVar.selectedResKind.push_back(kindIdx);
                break; // break from first loop
              }
            }
            // If there was no such a residue name, through error
            if (!found) {
              printf("Error: In Targeted-Swap move, residue name %s cannot be "
                     "swapped or was not found in PDB/PSF file!\n",
                     kindName.c_str());
              exit(EXIT_FAILURE);
            }
          }
        }
#if ENSEMBLE == GCMC
        // default value of chemPot is system value
        for (uint k = 0; k < molLookRef.GetNumKind(); k++) {
          tempVar.chemPot.push_back(molRef.kinds[k].chemPot);
        }
        // Set the chemical potential once we had a clear idea of our selected
        // reskind and overwrite the system value
        std::map<std::string, double>::const_iterator start =
            tsp.subVolumeChemPot.cp.begin();
        std::map<std::string, double>::const_iterator end =
            tsp.subVolumeChemPot.cp.end();
        while (start != end) {
          bool found = false;
          std::string resname = start->first;
          double cpot = start->second;
          for (int k = 0; k < (int)tempVar.selectedResKind.size(); ++k) {
            int kIndex = tempVar.selectedResKind[k];
            if (resname == molRef.kinds[kIndex].name) {
              found = true;
              tempVar.chemPot[kIndex] = cpot;
            }
          }
          // If there was no such a residue name, through error
          if (!found) {
            printf("Error: In Targeted-Swap move, residue name %s cannot be "
                   "swapped or was not found in PDB/PSF file!\n",
                   resname.c_str());
            exit(EXIT_FAILURE);
          }
          // make sure we use fugacity or chemical potential everywhere
          if (ffRef.isFugacity != tsp.subVolumeChemPot.isFugacity) {
            printf("Error: In Targeted-Swap move, chemPot and Fugacity cannot "
                   "be defined together!\n");
            exit(EXIT_FAILURE);
          }

          ++start;
        }
#endif
        targetSwapParam[tsp.selectedBox].push_back(tempVar);
        hasSubVolume[tsp.selectedBox] = true;
      }

      // Initialize the trial and acceptance
      for (uint b = 0; b < BOX_TOTAL; ++b) {
        if (hasSubVolume[b]) {
          trial[b].resize(targetSwapParam[b].size());
          accepted[b].resize(targetSwapParam[b].size());
          for (uint idx = 0; idx < targetSwapParam[b].size(); ++idx) {
            trial[b][idx].resize(molRef.GetKindsCount(), 0);
            accepted[b][idx].resize(molRef.GetKindsCount(), 0);
          }
        }
      }

      // Get the atom index for growing seed for each atomKind
      for (uint k = 0; k < molLookRef.GetNumKind(); k++) {
        growingAtomIndex.push_back(GetGrowingAtomIndex(k));
      }

      PrintTargetedSwapInfo();
    }
  }

  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0) {
    return mv::fail_state::NO_FAIL;
  }
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const ulong step);
  virtual void PrintAcceptKind();

private:
  double GetCoeff() const;
  uint GetBoxPairAndMol(const double subDraw, const double movPerc);
  uint PickMolInSubVolume();
  // Calculate the subvolume center, using defined atomList
  XYZ GetSubVolumeCenter(const uint &box, const uint &SubVIdx);
  // Check periodic boundary for subVolume
  void CheckSubVolumePBC(const uint &box);
  void PrintTargetedSwapInfo();
  // Track the acceptance for each subVolume
  void AcceptSubVolumeIdx(const uint &state, const uint &kind, const uint &box);
  // Returns the center node in molecule's graph
  uint GetGrowingAtomIndex(const uint k);
  // Finding the molecule of kind inside cavity and store the molecule Index
  bool FindMolInSubVolume(const uint box, const uint kind, const bool useGC);
  // find the molecule index that it's geometric center is within sub-volume
  template <bool pbcX, bool pbcY, bool pbcZ>
  bool SearchCavity_GC(std::vector<uint> &mol, const XYZ &center,
                       const XYZ &cavDim, const uint box, const uint kind);
  // find the molecule index that specific atom is within sub-volume
  template <bool pbcX, bool pbcY, bool pbcZ>
  bool SearchCavity_AC(std::vector<uint> &mol, const XYZ &center,
                       const XYZ &cavDim, const uint box, const uint kind,
                       const int atomIdx);
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  MoleculeLookup &molLookRef;
  Forcefield const &ffRef;

  std::vector<uint> growingAtomIndex;
  std::vector<TSwapParam> targetSwapParam[BOX_TOTAL];
  std::vector<uint> molIdxInSubVolume[BOX_TOTAL];
  // For move acceptance of each molecule subVolume, and kind
  std::vector<std::vector<uint>> trial[BOX_TOTAL], accepted[BOX_TOTAL];
  XYZ subVolCenter[BOX_TOTAL], subVolDim[BOX_TOTAL];
  double subVolume[BOX_TOTAL];
  bool hasSubVolume[BOX_TOTAL];
  int pickedSubV[BOX_TOTAL];
  // sub-Volume periodicity for each axis
  int pbcMode[BOX_TOTAL];
  bool rigidSwap;
  // sub-Volume periodicity for each axis
  std::vector<bool> xyzPBC[BOX_TOTAL]; // true if we have PBC in x, y, or z axis
};

// Need to fix it for GEMC
void TargetedSwap::PrintAcceptKind() {
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    if (hasSubVolume[b]) {
      for (uint idx = 0; idx < targetSwapParam[b].size(); ++idx) {
        for (uint k = 0; k < targetSwapParam[b][idx].selectedResKind.size();
             ++k) {
          int tsKind = targetSwapParam[b][idx].selectedResKind[k];
          int t = trial[b][idx][tsKind];
          int a = accepted[b][idx][tsKind];
          double percentAccept = (t ? 100 * (double)a / t : 0.0);
          char msg[257];
          sprintf(msg, "%% Accepted Targeted-Trans(%d)", idx);
          printf("%-30s %-5s ", msg, molRef.kinds[tsKind].name.c_str());
          printf("%10.5f %10.5f \n", percentAccept, percentAccept);
        }
      }
    }
  }
}

inline uint TargetedSwap::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_TARGETED_SWAP);
  overlap = false;
  // Pick the source, dest box. Pick the subVolume in each box (if exist)
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStart);

    // We need to copy the oldMol coordinate to newMol in order to be able to
    // insert rigidly in dest box. This means we need to unwrap it using source
    // box, and then wrap it using dest box
    XYZArray molA(pLen);
    coordCurrRef.CopyRange(molA, pStart, 0, pLen);
    boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndex));
    boxDimRef.WrapPBC(molA, destBox);
    // set coordinate of mole to newMol, later it will shift to center
    newMol.SetCoords(molA, 0);

    // the trial configuration has cavity but COM is not fix and no rotation
    // around backbone
    if (hasSubVolume[sourceBox]) { // deletion
      // Since the subVolume is orthogonal, no need to calculate and set the
      // cellBasis matrix. The default value for cellBasis matrix is I
      oldMol.SetSeed(subVolCenter[sourceBox], subVolDim[sourceBox], true, false,
                     false);
      oldMol.SetGrowingAtomIndex(growingAtomIndex[kindIndex]);
    }
    if (hasSubVolume[destBox]) { // insertion
      // Since the subVolume is orthogonal, no need to calculate and set the
      // cellBasis matrix. The default value for cellBasis matrix is I
      newMol.SetSeed(subVolCenter[destBox], subVolDim[destBox], true, false,
                     false);
      newMol.SetGrowingAtomIndex(growingAtomIndex[kindIndex]);
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_TARGETED_SWAP);
  return state;
}

inline uint TargetedSwap::GetBoxPairAndMol(const double subDraw,
                                           const double movPerc) {
  uint state = mv::fail_state::NO_FAIL;

  // 1. Pick a pair of box to determin source and destination box.
  prng.PickBoxPair(sourceBox, destBox, subDraw, movPerc);
  molIdxInSubVolume[sourceBox].clear();
  molIdxInSubVolume[destBox].clear();
  rigidSwap = true;
  // If we have a subvolume for source box, pick one of the subvolume randomely
  if (hasSubVolume[sourceBox]) {
    int SubVIdx_source = prng.randIntExc(targetSwapParam[sourceBox].size());
    // Set the picked subVolume parameter for source box
    pickedSubV[sourceBox] = SubVIdx_source;
    pbcMode[sourceBox] =
        targetSwapParam[sourceBox][SubVIdx_source].subVolumePBC;
    xyzPBC[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].xyzPBC;
    subVolDim[sourceBox] =
        targetSwapParam[sourceBox][SubVIdx_source].subVolumeDim;
    subVolume[sourceBox] = targetSwapParam[sourceBox][SubVIdx_source].subVolume;
    rigidSwap = targetSwapParam[sourceBox][SubVIdx_source].rigidSwap;
    // determin if we should calculate the center of subvolume or not
    if (targetSwapParam[sourceBox][SubVIdx_source].calcSubVolCenter) {
      subVolCenter[sourceBox] = GetSubVolumeCenter(sourceBox, SubVIdx_source);
    } else {
      subVolCenter[sourceBox] =
          targetSwapParam[sourceBox][SubVIdx_source].subVolumeCenter;
    }
    // Wrap the center of subVolume and make sure the it's dimension is less
    // than simulation box
    CheckSubVolumePBC(sourceBox);
  }

  // If we have a subvolume for dest box, pick one of the subvolume randomely
  if (hasSubVolume[destBox]) {
    int SubVIdx_dest = prng.randIntExc(targetSwapParam[destBox].size());
    // Set the picked subVolume parameter for dest box
    pickedSubV[destBox] = SubVIdx_dest;
    pbcMode[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumePBC;
    xyzPBC[destBox] = targetSwapParam[destBox][SubVIdx_dest].xyzPBC;
    subVolDim[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolumeDim;
    subVolume[destBox] = targetSwapParam[destBox][SubVIdx_dest].subVolume;
    rigidSwap = targetSwapParam[destBox][SubVIdx_dest].rigidSwap;
    // determin if we should calculate the center of subvolume or not
    if (targetSwapParam[destBox][SubVIdx_dest].calcSubVolCenter) {
      subVolCenter[destBox] = GetSubVolumeCenter(destBox, SubVIdx_dest);
    } else {
      subVolCenter[destBox] =
          targetSwapParam[destBox][SubVIdx_dest].subVolumeCenter;
    }
    // Wrap the center of subVolume and make sure the it's dimension is less
    // than simulation box
    CheckSubVolumePBC(destBox);
  }

  // 2. Pick a molecule kind
  // For the simulations, where all subVolumes are in one box we can do the
  // following. If there are subVolumes in each box, then we have to make sure
  // that targeted residue kind in all subvolumes are exactly identical.
  // Otherwise the results would
  //  would be biased.
  if (hasSubVolume[sourceBox]) {
    // Pick a molecule kind that specified for sourceBox subvolume
    uint i = prng.randIntExc(targetSwapParam[sourceBox][pickedSubV[sourceBox]]
                                 .selectedResKind.size());
    kindIndex =
        targetSwapParam[sourceBox][pickedSubV[sourceBox]].selectedResKind[i];
  } else if (hasSubVolume[destBox]) {
    // Pick a molecule kind that specified for destBox subvolume
    uint i = prng.randIntExc(
        targetSwapParam[destBox][pickedSubV[destBox]].selectedResKind.size());
    kindIndex =
        targetSwapParam[destBox][pickedSubV[destBox]].selectedResKind[i];
  } else {
    printf("Error: In Targeted-Swap move, no subVolume was defined for any "
           "Box!\n");
    exit(EXIT_FAILURE);
  }

  // 3. Pick a molecule Index for the picked molecule kind
  state = PickMolInSubVolume();

#if ENSEMBLE == GCMC
  if (state == mv::fail_state::NO_MOL_OF_KIND_IN_BOX && sourceBox == mv::BOX1 && hasSubVolume[sourceBox]) {
    std::cout << "Error: There are no molecules of kind "
              << molRef.kinds[kindIndex].name << " left in reservoir.\n";
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

inline XYZ TargetedSwap::GetSubVolumeCenter(const uint &box,
                                            const uint &subVIdx) {
  int i;
  TSwapParam tsp = targetSwapParam[box][subVIdx];
  int listSize = tsp.atomList.size();
  // use first atomIndex to unwrap others
  XYZ reference = coordCurrRef.Get(tsp.atomList[0]);
  XYZ center = reference;
  XYZ atomCoord;
  // start from next atom index
  for (i = 1; i < listSize; ++i) {
    if (molLookRef.IsAtomInBox(tsp.atomList[i], box)) {
      atomCoord = coordCurrRef.Get(tsp.atomList[i]);
      center += boxDimRef.UnwrapPBC(atomCoord, box, reference);
    } else {
      printf("Error: In Targeted-Swap move, atom index %d does not exist in "
             "box %d!\n",
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

inline void TargetedSwap::CheckSubVolumePBC(const uint &box) {
  XYZ boxAxis = boxDimRef.GetAxis(box);
  int subVolIdx = targetSwapParam[box][pickedSubV[box]].subVolumeIdx;
  bool shrunk = false;
  // Check the subVolume dimension
  if (subVolDim[box].x > boxAxis.x) {
    printf("Error: In Targeted-Swap move, box dimension in X axis shrunk "
           "bellow %g A for Subvolume %d!\n",
           subVolDim[box].x, subVolIdx);
    shrunk = true;
  }
  if (subVolDim[box].y > boxAxis.y) {
    printf("Error: In Targeted-Swap move, box dimension in Y axis shrunk "
           "bellow %g A for Subvolume %d!\n",
           subVolDim[box].y, subVolIdx);
    shrunk = true;
  }
  if (subVolDim[box].z > boxAxis.z) {
    printf("Error: In Targeted-Swap move, box dimension in Z axis shrunk "
           "bellow %g A for Subvolume %d!\n",
           subVolDim[box].z, subVolIdx);
    shrunk = true;
  }

  // Wrap the center of subVolume and check for PBC
  subVolCenter[box] = boxDimRef.WrapPBC(subVolCenter[box], box, xyzPBC[box][0],
                                        xyzPBC[box][1], xyzPBC[box][2]);
  // need to look at unslant center in order to compare with XYZ axis
  XYZ unSlanCenter = boxDimRef.TransformUnSlant(subVolCenter[box], box);
  // Check the center or subvolume
  if ((unSlanCenter.x >= boxAxis.x) || (unSlanCenter.x < 0)) {
    printf("Error: In Targeted-Swap move, subVolume center.x %g A for "
           "Subvolume %d is out of range of axis.x %g A!\n",
           subVolCenter[box].x, subVolIdx, boxAxis.x);
    printf("     : Define the subVolume center.x within the simulation box or "
           "turn on the PBC in x axis.\n");
    shrunk = true;
  }
  if ((unSlanCenter.y >= boxAxis.y) || (unSlanCenter.y < 0)) {
    printf("Error: In Targeted-Swap move, subVolume center.y %g A for "
           "Subvolume %d is out of range of axis.y %g A!\n",
           subVolCenter[box].y, subVolIdx, boxAxis.y);
    printf("     : Define the subVolume center.y within the simulation box or "
           "turn on the PBC in y axis.\n");
    shrunk = true;
  }
  if ((unSlanCenter.z >= boxAxis.z) || (unSlanCenter.z < 0)) {
    printf("Error: In Targeted-Swap move, subVolume center.z %g A for "
           "Subvolume %d is out of range of axis.z %g A!\n",
           subVolCenter[box].z, subVolIdx, boxAxis.z);
    printf("     : Define the subVolume center.z within the simulation box or "
           "turn on the PBC in z axis.\n");
    shrunk = true;
  }

  // check the subVolume dimensions with respect to simBox
  XYZ unSlantMaxDim = unSlanCenter + (subVolDim[box] * 0.5);
  XYZ unSlantMinDim = unSlanCenter - (subVolDim[box] * 0.5);
  XYZ slantMaxDim = subVolCenter[box] + (subVolDim[box] * 0.5);
  XYZ slantMinDim = subVolCenter[box] - (subVolDim[box] * 0.5);
  if (!xyzPBC[box][0]) { // no PBC in x axis
    // need to use >= because if atom is exactly on axis, it would wrap to 0.0!
    if ((unSlantMaxDim.x >= boxAxis.x) || (unSlantMinDim.x < 0)) {
      printf("Error: In Targeted-Swap move with no PBC in x axis, subVolume.x "
             "range [%g-%g] A for Subvolume %d is on the edge or out of range "
             "of axis.x %g A!\n",
             slantMinDim.x, slantMaxDim.x, subVolIdx, boxAxis.x);
      printf("     : Decrease the subVolume dimension.x or turn on the PBC in "
             "x axis.\n");
      shrunk = true;
    }
  }
  if (!xyzPBC[box][1]) { // no PBC in y axis
    // need to use >= because if atom is exactly on axis, it would wrap to 0.0!
    if ((unSlantMaxDim.y >= boxAxis.y) || (unSlantMinDim.y < 0)) {
      printf("Error: In Targeted-Swap move with no PBC in y axis, subVolume.y "
             "range [%g-%g] A for Subvolume %d is on the edge or out of range "
             "of axis.y %g A!\n",
             slantMinDim.y, slantMaxDim.y, subVolIdx, boxAxis.y);
      printf("     : Decrease the subVolume dimension.y or turn on the PBC in "
             "y axis.\n");
      shrunk = true;
    }
  }
  if (!xyzPBC[box][2]) { // no PBC in z axis
    // need to use >= because if atom is exactly on axis, it would wrap to 0.0!
    if ((unSlantMaxDim.z >= boxAxis.z) || (unSlantMinDim.z < 0)) {
      printf("Error: In Targeted-Swap move with no PBC in z axis, subVolume.z "
             "range [%g-%g] A for Subvolume %d is on the edge or out of range "
             "of axis.z %g A!\n",
             slantMinDim.z, slantMaxDim.z, subVolIdx, boxAxis.z);
      printf("     : Decrease the subVolume dimension.z or turn on the PBC in "
             "z axis.\n");
      shrunk = true;
    }
  }

  if (shrunk) {
    exit(EXIT_FAILURE);
  }
}

inline uint TargetedSwap::PickMolInSubVolume() {
  uint rejectState = mv::fail_state::NO_FAIL;
  if (hasSubVolume[sourceBox]) {
    // Search through all molecule in the subVolume and randomely pick
    // one of the picked kind
    bool foundMol = FindMolInSubVolume(sourceBox, kindIndex, rigidSwap);

    if (foundMol) {
      // The return vector, stores unique molecule index
      uint i = prng.randIntExc(molIdxInSubVolume[sourceBox].size());
      molIndex = molIdxInSubVolume[sourceBox][i];
    } else {
      // reject the move if no molecule was find
      rejectState = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    }
  } else { // means we have subVolume in destBox
    // Randomely pick one molecule of the picked kind from whole source box
    rejectState = prng.PickMolIndex(molIndex, kindIndex, sourceBox);
    // Just call the function to get number of molecule in cavity in destBox
    FindMolInSubVolume(destBox, kindIndex, rigidSwap);
  }

  return rejectState;
}

inline uint TargetedSwap::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_TARGETED_SWAP);
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  if (rigidSwap) {
    molRef.kinds[kindIndex].BuildIDOld(oldMol, molIndex);
    // Add bonded energy because we don't consider it in DCRotateCOM.cpp
    oldMol.AddEnergy(calcEnRef.MoleculeIntra(oldMol));

    molRef.kinds[kindIndex].BuildIDNew(newMol, molIndex);
    // Add bonded energy because we don't consider it in DCRotateCOM.cpp
    newMol.AddEnergy(calcEnRef.MoleculeIntra(newMol));
  } else {
    // insert or delete using CD-CBMC starting from growingatom Index
    molRef.kinds[kindIndex].BuildGrowInCav(oldMol, newMol, molIndex);
  }
  overlap = newMol.HasOverlap();
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_TARGETED_SWAP);
  return mv::fail_state::NO_FAIL;
}

inline void TargetedSwap::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_TARGETED_SWAP);
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
    // SwapDestRecip must be called first to backup the cosMol and sinMol
    recipGain.energy = calcEwald->SwapDestRecip(newMol, destBox, molIndex);
    recipLose.energy = calcEwald->SwapSourceRecip(oldMol, sourceBox, molIndex);
    // need to contribute the self and correction energy
    W_recip = exp(-1.0 * ffRef.beta *
                  (recipGain.energy + recipLose.energy + correct_new -
                   correct_old + self_new - self_old));
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_TARGETED_SWAP);
}

inline double TargetedSwap::GetCoeff() const {
#if ENSEMBLE == GEMC
  printf("Error: GOMC does not support Targeted-Swap for GEMC!\n");
  exit(EXIT_FAILURE);
  // Need to be implemented.
  // XXXXXXXXXXXXXXXXXXXXXX
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) { // Delete case
    double molNumber = double(molIdxInSubVolume[sourceBox].size());
    double invVolume = 1.0 / subVolume[sourceBox];
    double chemicalP =
        targetSwapParam[sourceBox][pickedSubV[sourceBox]].chemPot[kindIndex];
    /*
    if (molIdxInSubVolume[sourceBox].size() !=
    molLookRef.NumKindInBox(kindIndex, sourceBox)) { printf("Error: %d in
    subVolume vs %d in box!\n", molIdxInSubVolume[sourceBox].size(),
      molLookRef.NumKindInBox(kindIndex, sourceBox));
    } */
    if (ffRef.isFugacity) {
      return molNumber * invVolume / (BETA * chemicalP);
    } else {
      return molNumber * invVolume * exp(-BETA * chemicalP);
    }
  } else { // Insertion case
    double molNumber = double(molIdxInSubVolume[destBox].size());
    double volume = subVolume[destBox];
    double chemicalP =
        targetSwapParam[destBox][pickedSubV[destBox]].chemPot[kindIndex];

    /*
    if (molIdxInSubVolume[destBox].size() != molLookRef.NumKindInBox(kindIndex,
    destBox)) { printf("Error: %d in subVolume vs %d in box!\n",
    molIdxInSubVolume[destBox].size(), molLookRef.NumKindInBox(kindIndex,
    destBox));
    } */
    if (ffRef.isFugacity) {
      return volume / (molNumber + 1) * (BETA * chemicalP);
    } else {
      return volume / (molNumber + 1) * exp(BETA * chemicalP);
    }
  }
#endif
}

inline void TargetedSwap::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_TARGETED_SWAP);
  bool result;
  // If we didn't skip the move calculation
  if (rejectState == mv::fail_state::NO_FAIL) {
    double molTransCoeff = GetCoeff();
    double Wo = oldMol.GetWeight();
    double Wn = newMol.GetWeight();
    double Wrat = Wn / Wo * W_tc * W_recip;

    // safety to make sure move will be rejected in overlap case
    if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
      result = prng() < molTransCoeff * Wrat;
    } else
      result = false;

    if (result) {
      // Add tail corrections
      sysPotRef.boxEnergy[sourceBox].tailCorrection += tcLose.energy;
      sysPotRef.boxEnergy[destBox].tailCorrection += tcGain.energy;
      // Add rest of energy.
      sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
      sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
      //Add Reciprocal energy
      sysPotRef.boxEnergy[sourceBox].recip += recipLose.energy;
      sysPotRef.boxEnergy[destBox].recip += recipGain.energy;
      // Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
      sysPotRef.boxEnergy[destBox].correction += correct_new;
      // Add self energy
      sysPotRef.boxEnergy[sourceBox].self -= self_old;
      sysPotRef.boxEnergy[destBox].self += self_new;

      // Set coordinates, new COM; shift index to new box's list
      newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
      comCurrRef.SetNew(molIndex, destBox);
      molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
      cellList.AddMol(molIndex, destBox, coordCurrRef);

      // Zero out box energies to prevent small number
      // errors in double.
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

      // Retotal
      sysPotRef.Total();
      // Update the velocity
      velocity.UpdateMolVelocity(molIndex, destBox);
    } else {
      cellList.AddMol(molIndex, sourceBox, coordCurrRef);
      // when weight is 0, MolDestSwap() will not be executed, thus cos/sin
      // molRef will not be changed. Also since no memcpy, doing restore
      // results in memory overwrite
      if (newMol.GetWeight() != 0.0 && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }

  } else // we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::TARGETED_SWAP, result, destBox, kindIndex);
  AcceptSubVolumeIdx(result, kindIndex, destBox);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_TARGETED_SWAP);
}

uint TargetedSwap::GetGrowingAtomIndex(const uint k) {
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
    // Setup the node's degree
    uint bondCount = kind.bondList.count;
    for (int i = 0; i < (int)bondCount; ++i) {
      flw.AddEdge(kind.bondList.part1[i], kind.bondList.part2[i]);
    }

    // get the centric node
    int centric = flw.GetCentricNode();

    // If it was not successfull, it returns -1
    if (centric == -1) {
      printf("Error: In TargetedSwap move, no center atom was found for %s!\n",
             kind.name.c_str());
      exit(EXIT_FAILURE);
    }

    return centric;
  }
}

void TargetedSwap::PrintTargetedSwapInfo() {
  int i, k, b;
  for (b = 0; b < BOX_TOTAL; ++b) {
    for (i = 0; i < (int)targetSwapParam[b].size(); ++i) {
      TSwapParam tsp = targetSwapParam[b][i];
      printf("%-40s %d: \n",
             "Info: Targeted-Swap parameter for subVolume index",
             tsp.subVolumeIdx);
      printf("%-40s %d \n", "      SubVolume Box:", b);
      if (tsp.calcSubVolCenter) {
        printf("%-40s Using %zu defined atom indexes \n",
               "      Calculating subVolume center:", tsp.atomList.size());
        int max = *std::max_element(tsp.atomList.begin(), tsp.atomList.end());
        if (max >= (int)coordCurrRef.Count()) {
          printf("Error: Atom index %d is beyond total number of atoms (%d) in "
                 "the system!\n",
                 max, coordCurrRef.Count());
          printf("       Make sure to use 0 based atom index!\n");
          exit(EXIT_FAILURE);
        }
      } else {
        printf("%-40s %g %g %g \n",
               "      SubVolume center:", tsp.subVolumeCenter.x,
               tsp.subVolumeCenter.y, tsp.subVolumeCenter.z);
      }

      printf("%-40s %s%s%s \n",
             "      SubVolume PBC:", (tsp.xyzPBC[0] ? "X" : ""),
             (tsp.xyzPBC[1] ? "Y" : ""), (tsp.xyzPBC[2] ? "Z" : ""));
      printf("%-40s %g %g %g \n",
             "      SubVolume dimension:", tsp.subVolumeDim.x,
             tsp.subVolumeDim.y, tsp.subVolumeDim.z);
      printf("%-40s %s \n", "      SubVolume Swap type:",
             (tsp.rigidSwap ? "Rigid body" : "Flexible"));
      printf("%-40s ", "      Targeted residue kind:");
      for (k = 0; k < (int)tsp.selectedResKind.size(); ++k) {
        printf("%-5s ", molRef.kinds[tsp.selectedResKind[k]].name.c_str());
      }
      printf("\n");
      if (!tsp.rigidSwap) {
        for (k = 0; k < (int)tsp.selectedResKind.size(); ++k) {
          int kIdx = tsp.selectedResKind[k];
          printf("%-40s %s: (%d, %s) \n",
                 "      Starting atom (index, name) for",
                 molRef.kinds[kIdx].name.c_str(), growingAtomIndex[kIdx],
                 molRef.kinds[kIdx].atomNames[growingAtomIndex[kIdx]].c_str());
        }
      }
#if ENSEMBLE == GCMC
      printf("%-40s ",
             (ffRef.isFugacity ? "      Targeted fugacity (bar):"
                               : "      Targeted chemical potential (K):"));
      for (k = 0; k < (int)tsp.selectedResKind.size(); ++k) {
        int kIdx = tsp.selectedResKind[k];
        double value =
            tsp.chemPot[kIdx] /
            (ffRef.isFugacity ? unit::BAR_TO_K_MOLECULE_PER_A3 : 1.0);
        printf("(%-4s: %g) ", molRef.kinds[kIdx].name.c_str(), value);
      }
#endif
      printf("\n\n");
    }
  }
}

void TargetedSwap::AcceptSubVolumeIdx(const uint &state, const uint &kind,
                                      const uint &box) {
  if (hasSubVolume[box]) {
    uint subVidx = pickedSubV[box];
    ++trial[box][subVidx][kind];
    if (state) {
      ++accepted[box][subVidx][kind];
    }
  }
}

bool TargetedSwap::FindMolInSubVolume(const uint box, const uint kind,
                                      const bool useGC) {
  bool found = false;
  if (useGC) { // uses the geometric center to detect the molecule in cavity
    switch (pbcMode[box]) {
    case PBC_X:
      found = SearchCavity_GC<true, false, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_Y:
      found = SearchCavity_GC<false, true, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_Z:
      found = SearchCavity_GC<false, false, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_XY:
      found = SearchCavity_GC<true, true, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_XZ:
      found = SearchCavity_GC<true, false, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_YZ:
      found = SearchCavity_GC<false, true, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;
    case PBC_XYZ:
      found = SearchCavity_GC<true, true, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind);
      break;

    default:
      printf("Error: Unknown PBC mode %d in targetedSwap move!\n",
             pbcMode[box]);
      exit(EXIT_FAILURE);
      break;
    }

  } else { // uses the specified atom index to detect the molecule in cavity
    switch (pbcMode[box]) {
    case PBC_X:
      found = SearchCavity_AC<true, false, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_Y:
      found = SearchCavity_AC<false, true, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_Z:
      found = SearchCavity_AC<false, false, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_XY:
      found = SearchCavity_AC<true, true, false>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_XZ:
      found = SearchCavity_AC<true, false, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_YZ:
      found = SearchCavity_AC<false, true, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;
    case PBC_XYZ:
      found = SearchCavity_AC<true, true, true>(
          molIdxInSubVolume[box], subVolCenter[box], subVolDim[box], box, kind,
          growingAtomIndex[kind]);
      break;

    default:
      printf("Error: Unknown PBC mode %d in targetedSwap move!\n",
             pbcMode[box]);
      exit(EXIT_FAILURE);
      break;
    }
  }

  return found;
}

template <bool pbcX, bool pbcY, bool pbcZ>
bool TargetedSwap::SearchCavity_GC(std::vector<uint> &mol, const XYZ &center,
                                   const XYZ &cavDim, const uint box,
                                   const uint kind) {
  int *particleMol = molLookRef.molIndex;
  uint molKind, molIndex;
  XYZ halfDim = cavDim * 0.5;
  double maxLength = halfDim.Max();
  XYZ halfDimSq = halfDim * halfDim;
  XYZ dist, distSq; // distance from center to geometric center of molecule

  if (maxLength <= boxDimRef.rCut[box]) {
    CellList::Neighbors n = cellList.EnumerateLocal(center, box);
    while (!n.Done()) {
      molIndex = particleMol[*n];
      molKind = molRef.GetMolKind(molIndex);
      // Check the kind to before calculating distance to save time
      // if molecule can be transfer between boxes and it's the right kind
      if (!molLookRef.IsNoSwap(molIndex) && (molKind == kind)) {
        dist = comCurrRef.Get(molIndex) - center;
        // Apply pbc on selectice axis
        if (pbcX) {
          dist = boxDimRef.MinImage_X(dist, box);
        }
        if (pbcY) {
          dist = boxDimRef.MinImage_Y(dist, box);
        }
        if (pbcZ) {
          dist = boxDimRef.MinImage_Z(dist, box);
        }
        distSq = dist * dist;
        if (distSq.x <= halfDimSq.x && distSq.y <= halfDimSq.y &&
            distSq.z <= halfDimSq.z) {
          mol.push_back(molIndex);
        }
      }
      n.Next();
    }

    // Find a unique molecule index
    std::vector<uint>::iterator ip;
    std::sort(mol.begin(), mol.end());
    ip = std::unique(mol.begin(), mol.end());
    mol.resize(std::distance(mol.begin(), ip));
  } else {
    MoleculeLookup::box_iterator n = molLookRef.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookRef.BoxEnd(box);
    // uint start, length, p;
    while (n != end) {
      molIndex = *n;
      molKind = molRef.GetMolKind(molIndex);
      // Check the kind to before calculating distance to save time
      // if molecule can be transfer between boxes and it's the right kind
      if (!molLookRef.IsNoSwap(molIndex) && (molKind == kind)) {
        dist = comCurrRef.Get(molIndex) - center;
        // Apply pbc on selectice axis
        if (pbcX) {
          dist = boxDimRef.MinImage_X(dist, box);
        }
        if (pbcY) {
          dist = boxDimRef.MinImage_Y(dist, box);
        }
        if (pbcZ) {
          dist = boxDimRef.MinImage_Z(dist, box);
        }

        distSq = dist * dist;
        if (distSq.x <= halfDimSq.x && distSq.y <= halfDimSq.y &&
            distSq.z <= halfDimSq.z) {
          mol.push_back(molIndex);
        }
      }
      n++;
    }
    // No need to find the unique molIndex, since we loop through all molecules
    // and and not atoms
  }

  // returns true if there is a molecule of kind in cavity
  return mol.size();
}

template <bool pbcX, bool pbcY, bool pbcZ>
bool TargetedSwap::SearchCavity_AC(std::vector<uint> &mol, const XYZ &center,
                                   const XYZ &cavDim, const uint box,
                                   const uint kind, const int atomIdx) {
  int *particleMol = molLookRef.molIndex;
  int *particleIndex = molLookRef.atomIndex;
  uint molKind, molIndex;
  int aIdx;
  XYZ halfDim = cavDim * 0.5;
  double maxLength = halfDim.Max();
  XYZ halfDimSq = halfDim * halfDim;
  XYZ dist, distSq; // distance from center to atoms

  if (maxLength <= boxDimRef.rCut[box]) {
    CellList::Neighbors n = cellList.EnumerateLocal(center, box);
    while (!n.Done()) {
      aIdx = particleIndex[*n];
      molIndex = particleMol[*n];
      molKind = molRef.GetMolKind(molIndex);
      // Check the kind to before calculating distance to save time
      // if molecule can be transfer between boxes and it's the right kind
      if (aIdx == atomIdx && (molKind == kind) &&
          !molLookRef.IsNoSwap(molIndex)) {
        dist = coordCurrRef.Get(*n) - center;
        // Apply pbc on selectice axis
        if (pbcX) {
          dist = boxDimRef.MinImage_X(dist, box);
        }
        if (pbcY) {
          dist = boxDimRef.MinImage_Y(dist, box);
        }
        if (pbcZ) {
          dist = boxDimRef.MinImage_Z(dist, box);
        }

        distSq = dist * dist;
        if (distSq.x <= halfDimSq.x && distSq.y <= halfDimSq.y &&
            distSq.z <= halfDimSq.z) {
          mol.push_back(molIndex);
        }
      }
      n.Next();
    }

    // There should be only one atom with atomIdx, so no need to
    // find the unique molIndex
  } else {
    MoleculeLookup::box_iterator n = molLookRef.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookRef.BoxEnd(box);
    int start, length, p;
    while (n != end) {
      molIndex = *n;
      molKind = molRef.GetMolKind(molIndex);
      // Check the kind to before calculating distance to save time
      // if molecule can be transfer between boxes and it's the right kind
      if (!molLookRef.IsNoSwap(molIndex) && (molKind == kind)) {
        // get atom start index and length
        length = molRef.GetKind(molIndex).NumAtoms();
        start = molRef.MolStart(molIndex);
        // loop through all atom and calculate the distance between center and
        // atomIdx coordinate
        for (p = 0; p < length; ++p) {
          if (p == atomIdx) {
            dist = coordCurrRef.Get(start + p) - center;
            // Apply pbc on selectice axis
            if (pbcX) {
              dist = boxDimRef.MinImage_X(dist, box);
            }
            if (pbcY) {
              dist = boxDimRef.MinImage_Y(dist, box);
            }
            if (pbcZ) {
              dist = boxDimRef.MinImage_Z(dist, box);
            }

            distSq = dist * dist;
            if (distSq.x <= halfDimSq.x && distSq.y <= halfDimSq.y &&
                distSq.z <= halfDimSq.z) {
              mol.push_back(molIndex);
            }
            break;
          }
        }
      }
      n++;
    }
    // No need to find the unique molIndex, since we loop through all molecules
    // and and not atoms
  }

  // returns true if there is a molecule of kind in cavity
  return mol.size();
}

#endif

#endif
