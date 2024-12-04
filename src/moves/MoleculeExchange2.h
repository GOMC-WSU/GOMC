/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLECULEEXCHANGE2_H
#define MOLECULEEXCHANGE2_H

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
#include "GeomLib.h"
#include "MoleculeExchange1.h"
#include "TrialMol.h"

using std::vector;
using namespace geom;

// MEMC-2 Move:
//
// Swapping one Large molecule with one or more small molecules in dense phase
// and vice versa.
// Sub-Volume location and orientation is based on the COM and backbone of the
// the small molecule.

class MoleculeExchange2 : public MoleculeExchange1 {
public:
  MoleculeExchange2(System &sys, StaticVals const &statV)
      : MoleculeExchange1(sys, statV) {
    if (enableID) {
      SetMEMC(statV);
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

protected:
  virtual void AdjustExRatio();
  virtual void SetMEMC(StaticVals const &statV);
  virtual void SetExchangeData();
  virtual uint PickMolInCav();
  virtual uint ReplaceMolecule();
  virtual double GetCoeff() const;

  int smallBB[2];
  // To store total sets of exchange pairs
  vector<vector<uint>> smallBBVec;
};

inline void MoleculeExchange2::SetMEMC(StaticVals const &statV) {
  for (uint t = 0; t < exchangeRatioVec.size(); t++) {
    smallBB[0] = smallBB[1] = -1;
    for (uint i = 0; i < molRef.kinds[kindSVec[t]].NumAtoms(); i++) {
      if (molRef.kinds[kindSVec[t]].atomNames[i] ==
          statV.memcVal.smallBBAtom1[t]) {
        smallBB[0] = i;
      }
      if (molRef.kinds[kindSVec[t]].atomNames[i] ==
          statV.memcVal.smallBBAtom2[t]) {
        smallBB[1] = i;
      }
    }

    for (uint i = 0; i < 2; i++) {
      if (smallBB[i] == -1) {
        printf("Error: In ME-2 move, atom name %s or %s was not found in %s "
               "residue.\n",
               statV.memcVal.smallBBAtom1[t].c_str(),
               statV.memcVal.smallBBAtom2[t].c_str(),
               statV.memcVal.smallKind[t].c_str());
        exit(EXIT_FAILURE);
      }
    }

    if (molRef.kinds[kindSVec[t]].NumAtoms() > 1) {
      if (smallBB[0] == smallBB[1]) {
        printf("Error: In ME-2 move, atom names in small molecule backbone "
               "cannot be same!\n");
        exit(EXIT_FAILURE);
      }
    }
    vector<uint> temp(smallBB, smallBB + 2);
    smallBBVec.push_back(temp);
  }
}

inline void MoleculeExchange2::AdjustExRatio() {
  if (((counter + 1) % perAdjust) == 0) {
    int exMax = ceil((float)molInCavCount / (float)perAdjust);
    int exMin = 1;
    uint index = kindS + kindL * molRef.GetKindsCount();
    double currAccept = (double)(accepted[sourceBox][index]) /
                        (double)(trial[sourceBox][index]);
    if (std::abs(currAccept - lastAccept) >= 0.05 * currAccept) {
      if (currAccept >= lastAccept) {
        exchangeRatio += exDiff;
      } else {
        exDiff *= -1;
        exchangeRatio += exDiff;
      }
      lastAccept = currAccept;
      if (exchangeRatio < exMin)
        exchangeRatio = exMin;
      if (exchangeRatio > exMax)
        exchangeRatio = exMax;
    }
    molInCavCount = 0;
    counter = 0;
  }
}

inline void MoleculeExchange2::SetExchangeData() {
  uint exType = prng.randIntExc(exchangeRatioVec.size());
  kindS = kindSVec[exType];
  kindL = kindLVec[exType];
  exchangeRatio = exchangeRatioVec[exType];
  largeBB[0] = largeBBVec[exType][0];
  largeBB[1] = largeBBVec[exType][1];
  smallBB[0] = smallBBVec[exType][0];
  smallBB[1] = smallBBVec[exType][1];
}

inline uint MoleculeExchange2::PickMolInCav() {
  uint state = mv::fail_state::NO_FAIL;
  // pick a random small kind in dense phase and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if (state == mv::fail_state::NO_FAIL) {
    center = comCurrRef.Get(pickedS);
    // If we want to orient the cavity with backbone of picked small mol
    if (molRef.NumAtoms(kindS) == 1) {
      XYZ cavVect = prng.RandomUnitVect();
      SetBasis(cavA, cavVect);
    } else {
      uint start = molRef.MolStart(pickedS) + smallBB[0];
      uint end = molRef.MolStart(pickedS) + smallBB[1];
      SetBasis(cavA, boxDimRef.MinImage(coordCurrRef.Difference(start, end),
                                        sourceBox));
    }
    // Calculate inverse matrix for cav here Inv = transpose
    TransposeMatrix(invCavA, cavA);

    // Find the molecule kind 0 in the cavity
    if ((exchangeRatio == 1) ||
        calcEnRef.FindMolInCavity(molInCav, center, cavity, invCavA, sourceBox,
                                  kindS, exchangeRatio)) {
      // Find the exchangeRatio number of molecules kind 0 in cavity
      numInCavA = exchangeRatio;
      // add the random picked small molecule to the list.
      molIndexA.push_back(pickedS);
      kindIndexA.push_back(pickedKS);
      if (exchangeRatio == 1) {
        totMolInCav = 1;
      } else {
        totMolInCav = molInCav[kindS].size();
        // delete the picked small molecule from list
        for (uint s = 0; s < totMolInCav; s++) {
          if (pickedS == molInCav[kindS][s])
            molInCav[kindS].erase(molInCav[kindS].begin() + s);
        }
      }

      for (uint n = 1; n < numInCavA; n++) {
        // pick random exchangeRatio number of kindS in cavity
        uint picked = prng.randIntExc(molInCav[kindS].size());
        molIndexA.push_back(molInCav[kindS][picked]);
        kindIndexA.push_back(molRef.GetMolKind(molIndexA[n]));
        molInCav[kindS].erase(molInCav[kindS].begin() + picked);
      }

      // pick a molecule from Large kind in destBox
      numInCavB = 1;
      state = prng.PickMol(kindL, kindIndexB, molIndexB, numInCavB, destBox);
    } else {
      // reject the move
      state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    }

    molInCavCount += totMolInCav;
    counter++;
  }

  return state;
}

inline uint MoleculeExchange2::ReplaceMolecule() {
  uint state = mv::fail_state::NO_FAIL;
  numInCavA = 1;
  numInCavB = exchangeRatio;
  // pick a random molecule of Large kind in dens box
  state = prng.PickMol(kindL, kindIndexA, molIndexA, numInCavA, sourceBox);

  if (state == mv::fail_state::NO_FAIL) {
    // Set the V1 to the backbone of the large molecule
    if (molRef.NumAtoms(kindL) == 1) {
      XYZ cavVect = prng.RandomUnitVect();
      SetBasis(cavA, cavVect);
    } else {
      uint start = molRef.MolStart(molIndexA[0]) + largeBB[0];
      uint end = molRef.MolStart(molIndexA[0]) + largeBB[1];
      SetBasis(cavA, boxDimRef.MinImage(coordCurrRef.Difference(start, end),
                                        sourceBox));
    }
    // Calculate inverse matrix for cav. Here Inv = Transpose
    TransposeMatrix(invCavA, cavA);
    // Use to shift to the COM of new molecule
    center = comCurrRef.Get(molIndexA[0]);
    if (exchangeRatio == 1) {
      totMolInCav = 0;
    } else {
      // find how many of KindS exist in this center
      calcEnRef.FindMolInCavity(molInCav, center, cavity, invCavA, sourceBox,
                                kindS, exchangeRatio);
      totMolInCav = molInCav[kindS].size();
    }
    // pick exchangeRatio number of Small molecule from dest box
    state = prng.PickMol(kindS, kindIndexB, molIndexB, numInCavB, destBox);
  }
  return state;
}

inline uint MoleculeExchange2::Prep(const double subDraw,
                                    const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MEMC);
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    numTypeASource =
        (double)(molLookRef.NumKindInBoxSwappable(kindIndexA[0], sourceBox));
    numTypeADest =
        (double)(molLookRef.NumKindInBoxSwappable(kindIndexA[0], destBox));
    numTypeBSource =
        (double)(molLookRef.NumKindInBoxSwappable(kindIndexB[0], sourceBox));
    numTypeBDest =
        (double)(molLookRef.NumKindInBoxSwappable(kindIndexB[0], destBox));
    // transferring type A from source to dest
    for (uint n = 0; n < numInCavA; n++) {
      newMolA.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef, destBox));
      oldMolA.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef, sourceBox));
    }

    for (uint n = 0; n < numInCavB; n++) {
      // transferring type B from dest to source
      newMolB.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef, sourceBox));
      oldMolB.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef, destBox));
    }

    // set the old coordinate and new after proper wrap & unwrap
    for (uint n = 0; n < numInCavA; n++) {
      XYZArray molA(pLenA[n]);
      coordCurrRef.CopyRange(molA, pStartA[n], 0, pLenA[n]);
      boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA[n]));
      boxDimRef.WrapPBC(molA, destBox);
      oldMolA[n].SetCoords(coordCurrRef, pStartA[n]);
      // set coordinate of moleA to newMolA, later it will shift to center
      newMolA[n].SetCoords(molA, 0);
      // copy cavA matrix to slant the old trial of molA
      oldMolA[n].SetCavMatrix(cavA);
    }

    for (uint n = 0; n < numInCavB; n++) {
      XYZArray molB(pLenB[n]);
      coordCurrRef.CopyRange(molB, pStartB[n], 0, pLenB[n]);
      boxDimRef.UnwrapPBC(molB, destBox, comCurrRef.Get(molIndexB[n]));
      boxDimRef.WrapPBC(molB, sourceBox);
      oldMolB[n].SetCoords(coordCurrRef, pStartB[n]);
      // set coordinate of moleB to newMolB, later it will shift
      newMolB[n].SetCoords(molB, 0);
      // copy cavA matrix to slant the new trial of molB
      newMolB[n].SetCavMatrix(cavA);
    }

    for (uint n = 0; n < numInCavB; n++) {
      // SetSeed(has cavity, COM is fixed, rotate around Backbone)
      if (insertL) {
        // Inserting Lmol from destBox to the center of cavity in sourceBox
        newMolB[n].SetSeed(center, cavity, true, true, true);
        // Set the Backbone of large molecule to be inserted
        newMolB[n].SetBackBone(largeBB);
        // perform rotational trial move in destBox for L oldMol
        oldMolB[n].SetSeed(false, false, false);
      } else {
        if (n == 0) {
          // Inserting Small from destBox to the center of cavity in sourceBox
          newMolB[n].SetSeed(center, cavity, true, true, true);
          // Set the Backbone of small molecule to be inserted
          newMolB[n].SetBackBone(smallBB);
        } else {
          // Inserting S mol from destBox to the cavity in sourceBox
          newMolB[n].SetSeed(center, cavity, true, false, false);
        }
        // perform trial move in destBox for S oldMol
        oldMolB[n].SetSeed(false, false, false);
      }
    }

    for (uint n = 0; n < numInCavA; n++) {
      if (insertL) {
        // Inserting S mol from sourceBox to destBox
        newMolA[n].SetSeed(false, false, false);
        if (n == 0) {
          // perform trial move in cavity with fix COM for S oldMol
          oldMolA[n].SetSeed(center, cavity, true, true, true);
          // Set the Backbone of small molecule to be deleted
          oldMolA[n].SetBackBone(smallBB);
        } else {
          // perform trial move in cavity in sourceBox for S oldMol
          oldMolA[n].SetSeed(center, cavity, true, false, false);
        }
      } else {
        // Inserting L mol from sourceBox to destBox
        newMolA[n].SetSeed(false, false, false);
        // perform rotational trial move on COM for L oldMol
        oldMolA[n].SetSeed(center, cavity, true, true, true);
        // Set the Backbone of large molecule to be deleted
        oldMolA[n].SetBackBone(largeBB);
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MEMC);
  return state;
}

inline uint MoleculeExchange2::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MEMC);
  CalcTc();

  // Calc Old energy and delete A from source
  if (insertL) {
    // Remove the fixed COM small mol at the end because we insert it at first
    for (uint n = numInCavA; n > 0; n--) {
      cellList.RemoveMol(molIndexA[n - 1], sourceBox, coordCurrRef);
      molRef.kinds[kindIndexA[n - 1]].BuildIDOld(oldMolA[n - 1],
                                                 molIndexA[n - 1]);
      // Add bonded energy because we don't consider it in DCRotate.cpp
      oldMolA[n - 1].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n - 1]));
    }
  } else {
    for (uint n = 0; n < numInCavA; n++) {
      cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
      molRef.kinds[kindIndexA[n]].BuildIDOld(oldMolA[n], molIndexA[n]);
      // Add bonded energy because we don't consider it in DCRotate.cpp
      oldMolA[n].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n]));
    }
  }

  // Calc old energy and delete B from destBox
  for (uint n = 0; n < numInCavB; n++) {
    cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
    molRef.kinds[kindIndexB[n]].BuildIDOld(oldMolB[n], molIndexB[n]);
    // Add bonded energy because we don't consider it in DCRotate.cpp
    oldMolB[n].AddEnergy(calcEnRef.MoleculeIntra(oldMolB[n]));
  }

  // Insert A to destBox
  for (uint n = 0; n < numInCavA; n++) {
    molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
    ShiftMol(true, n, sourceBox, destBox);
    cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
    // Add bonded energy because we don't consider it in DCRotate.cpp
    newMolA[n].AddEnergy(calcEnRef.MoleculeIntra(newMolA[n]));
    overlap |= newMolA[n].HasOverlap();
  }

  // Insert B in sourceBox
  for (uint n = 0; n < numInCavB; n++) {
    molRef.kinds[kindIndexB[n]].BuildIDNew(newMolB[n], molIndexB[n]);
    ShiftMol(false, n, destBox, sourceBox);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
    // Add bonded energy because we don't consider it in DCRotate.cpp
    newMolB[n].AddEnergy(calcEnRef.MoleculeIntra(newMolB[n]));
    overlap |= newMolB[n].HasOverlap();
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MEMC);
  return mv::fail_state::NO_FAIL;
}

inline void MoleculeExchange2::CalcEn() { MoleculeExchange1::CalcEn(); }

inline double MoleculeExchange2::GetCoeff() const {
#if ENSEMBLE == GEMC
  double volDest = boxDimRef.volume[destBox];
  if (insertL) {
    // kindA is the small molecule
    double ratioF = num::Factorial(totMolInCav - 1) /
                    (num::Factorial(totMolInCav - exchangeRatio) *
                     num::Factorial(numTypeADest, exchangeRatio));
    double ratioV = pow(volDest / volCav, exchangeRatio - 1);
    double ratioM = numTypeASource * numTypeBDest / (numTypeBSource + 1.0);
    return ratioF * ratioV * ratioM;
  } else {
    // kindA is the big molecule
    double ratioF =
        num::Factorial(totMolInCav) *
        num::Factorial(numTypeBDest - exchangeRatio, exchangeRatio) /
        num::Factorial(totMolInCav + exchangeRatio - 1);
    double ratioV = pow(volCav / volDest, exchangeRatio - 1);
    double ratioM = numTypeASource /
                    ((numTypeADest + 1.0) * (numTypeBSource + exchangeRatio));
    return ratioF * ratioV * ratioM;
  }
#elif ENSEMBLE == GCMC
  if (ffRef.isFugacity) {
    double delA = molRef.kinds[kindIndexA[0]].chemPot * numInCavA;
    double insB = molRef.kinds[kindIndexB[0]].chemPot * numInCavB;
    if (insertL) {
      // Insert Large molecule
      double ratioF = num::Factorial(totMolInCav - 1) /
                      num::Factorial(totMolInCav - exchangeRatio);
      double ratioM = numTypeASource / (numTypeBSource + 1.0);
      return (insB / delA) * ratioF * ratioM / pow(volCav, exchangeRatio - 1);
    } else {
      // Delete Large Molecule
      double ratioF = num::Factorial(totMolInCav) /
                      num::Factorial(totMolInCav + exchangeRatio - 1);
      double ratioM = numTypeASource / (numTypeBSource + exchangeRatio);
      return (insB / delA) * ratioF * ratioM * pow(volCav, exchangeRatio - 1);
    }
  } else {
    double delA = (-BETA * molRef.kinds[kindIndexA[0]].chemPot * numInCavA);
    double insB = (BETA * molRef.kinds[kindIndexB[0]].chemPot * numInCavB);
    if (insertL) {
      // Insert Large molecule
      double ratioF = num::Factorial(totMolInCav - 1) /
                      num::Factorial(totMolInCav - exchangeRatio);
      double ratioM = numTypeASource / (numTypeBSource + 1.0);
      return exp(delA + insB) * ratioF * ratioM /
             pow(volCav, exchangeRatio - 1);
    } else {
      // Delete Large Molecule
      double ratioF = num::Factorial(totMolInCav) /
                      num::Factorial(totMolInCav + exchangeRatio - 1);
      double ratioM = numTypeASource / (numTypeBSource + exchangeRatio);
      return exp(delA + insB) * ratioF * ratioM *
             pow(volCav, exchangeRatio - 1);
    }
  }
#endif
}

inline void MoleculeExchange2::Accept(const uint rejectState,
                                      const ulong step) {
  MoleculeExchange1::Accept(rejectState, step);
}

#endif

#endif
