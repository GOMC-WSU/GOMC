/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLECULEEXCHANGE3_H
#define MOLECULEEXCHANGE3_H

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC

#include <cmath>

#include "GeomLib.h"
#include "MoleculeExchange1.h"
#include "TrialMol.h"

using std::vector;
using namespace geom;

// MEMC-3 Move:
//
// Swapping one Large molecule with one or more small molecules in dense phase
// and vice versa.
// Sub-Volume location at the COM of the small molecule with random orientation.
// Large molecule insertion and deletion is performed by CD-CBMC.

class MoleculeExchange3 : public MoleculeExchange1 {
public:
  MoleculeExchange3(System &sys, StaticVals const &statV)
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
};

inline void MoleculeExchange3::SetMEMC(StaticVals const &statV) {
  for (uint t = 0; t < exchangeRatioVec.size(); t++) {
    if (largeBBVec[t][0] != largeBBVec[t][1]) {
      printf("Error: In ME-3 move, two atoms with same name should be used as "
             "backbone.\n");
      printf("Atom names in backbone was set to %s or %s in %s residue.\n",
             statV.memcVal.largeBBAtom1[t].c_str(),
             statV.memcVal.largeBBAtom2[t].c_str(),
             statV.memcVal.largeKind[t].c_str());
      exit(EXIT_FAILURE);
    }
  }
}

inline void MoleculeExchange3::AdjustExRatio() {
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

inline void MoleculeExchange3::SetExchangeData() {
  MoleculeExchange1::SetExchangeData();
}

inline uint MoleculeExchange3::PickMolInCav() {
  uint state = mv::fail_state::NO_FAIL;
  // pick a random small kind in dense phase and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if (state == mv::fail_state::NO_FAIL) {
    center = comCurrRef.Get(pickedS);
    // Pick random vector and find two vectors that are perpendicular to V1
    XYZ cavVect = prng.RandomUnitVect();
    SetBasis(cavA, cavVect);
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

inline uint MoleculeExchange3::ReplaceMolecule() {
  uint state = mv::fail_state::NO_FAIL;
  numInCavA = 1;
  numInCavB = exchangeRatio;
  // pick a random molecule of Large kind in dens box
  state = prng.PickMol(kindL, kindIndexA, molIndexA, numInCavA, sourceBox);

  if (state == mv::fail_state::NO_FAIL) {
    // Set V1 to a random vector and calculate two vector perpendicular to V1
    XYZ cavVect = prng.RandomUnitVect();
    SetBasis(cavA, cavVect);
    // Calculate inverse matrix for cav. Here Inv = Transpose
    TransposeMatrix(invCavA, cavA);
    // use the predefine atom in kindL as the center
    uint start = molRef.MolStart(molIndexA[0]) + largeBB[0];
    center = coordCurrRef.Get(start);
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

inline uint MoleculeExchange3::Prep(const double subDraw,
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
        // Set the atom of large molecule to be inserted in COM of cavity
        newMolB[n].SetBackBone(largeBB);
        // perform rotational trial move in destBox for L oldMol
        oldMolB[n].SetSeed(false, false, false);
      } else {
        if (n == 0) {
          // Inserting Smol from destBox to the center of cavity in sourceBox
          newMolB[n].SetSeed(center, cavity, true, true, false);
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
          oldMolA[n].SetSeed(center, cavity, true, true, false);
        } else {
          // perform trial move in cavity in sourceBox for S oldMol
          oldMolA[n].SetSeed(center, cavity, true, false, false);
        }
      } else {
        // Inserting L mol from sourceBox to destBox
        newMolA[n].SetSeed(false, false, false);
        // perform rotational trial move on COM for L oldMol
        oldMolA[n].SetSeed(center, cavity, true, true, true);
        // Set the atom of the large molecule to be inserted in COM
        oldMolA[n].SetBackBone(largeBB);
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MEMC);
  return state;
}

inline uint MoleculeExchange3::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MEMC);
  CalcTc();

  // Deleting A, B from their box
  if (insertL) {
    // Remove the fixed COM small mol at the end because we insert it at first
    for (uint n = numInCavA; n > 0; n--) {
      cellList.RemoveMol(molIndexA[n - 1], sourceBox, coordCurrRef);
      molRef.kinds[kindIndexA[n - 1]].BuildIDOld(oldMolA[n - 1],
                                                 molIndexA[n - 1]);
      // Add bonded energy because we don't consider it in DCRotate.cpp
      oldMolA[n - 1].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n - 1]));
    }
    // Calc old energy and delete Large kind from dest box
    for (uint n = 0; n < numInCavB; n++) {
      cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
      molRef.kinds[kindIndexB[n]].BuildOld(oldMolB[n], molIndexB[n]);
    }
  } else {
    // Calc old energy and delete Large kind from source box
    for (uint n = 0; n < numInCavA; n++) {
      cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
      molRef.kinds[kindIndexA[n]].BuildGrowOld(oldMolA[n], molIndexA[n]);
    }
    // Calc old energy and delete Small kind from dest box
    for (uint n = 0; n < numInCavB; n++) {
      cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
      molRef.kinds[kindIndexB[n]].BuildIDOld(oldMolB[n], molIndexB[n]);
      oldMolB[n].AddEnergy(calcEnRef.MoleculeIntra(oldMolB[n]));
    }
  }

  // Inserting A, B to new box
  if (insertL) {
    // Insert Small kind to destBox
    for (uint n = 0; n < numInCavA; n++) {
      molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
      ShiftMol(true, n, sourceBox, destBox);
      cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
      newMolA[n].AddEnergy(calcEnRef.MoleculeIntra(newMolA[n]));
      overlap |= newMolA[n].HasOverlap();
    }
    // Insert Large kind to sourceBox
    for (uint n = 0; n < numInCavB; n++) {
      molRef.kinds[kindIndexB[n]].BuildGrowNew(newMolB[n], molIndexB[n]);
      ShiftMol(false, n, destBox, sourceBox);
      cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
      overlap |= newMolB[n].HasOverlap();
    }
  } else {
    // Insert Large kind to destBox
    for (uint n = 0; n < numInCavA; n++) {
      molRef.kinds[kindIndexA[n]].BuildNew(newMolA[n], molIndexA[n]);
      ShiftMol(true, n, sourceBox, destBox);
      cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
      overlap |= newMolA[n].HasOverlap();
    }
    // Insert Small kind to sourceBox
    for (uint n = 0; n < numInCavB; n++) {
      molRef.kinds[kindIndexB[n]].BuildIDNew(newMolB[n], molIndexB[n]);
      ShiftMol(false, n, destBox, sourceBox);
      cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
      // Add bonded energy because we don't consider it in DCRotate.cpp
      newMolB[n].AddEnergy(calcEnRef.MoleculeIntra(newMolB[n]));
      overlap |= newMolB[n].HasOverlap();
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MEMC);
  return mv::fail_state::NO_FAIL;
}

inline void MoleculeExchange3::CalcEn() { MoleculeExchange1::CalcEn(); }

inline double MoleculeExchange3::GetCoeff() const {
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

inline void MoleculeExchange3::Accept(const uint rejectState,
                                      const ulong step) {
  MoleculeExchange1::Accept(rejectState, step);
}

#endif

#endif
