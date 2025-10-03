/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef INTRAMOLECULEEXCHANGE3_H
#define INTRAMOLECULEEXCHANGE3_H

#include "GeomLib.h"
#include "IntraMoleculeExchange1.h"
#include "TrialMol.h"

using std::vector;
using namespace geom;

// Intra Molecule Exchange Move:
// KindA is small kind. KindB is large kind
// Orientation center of cavA is on COM of kindS, random orientation.
// Orientation center of cavB is on COM of kindL, aligned with kindL backbone.
// Delete the exchangeRatio kindS from cavA, and 1 kindL from cavB.
// Insert the exchangeRatio kindS to cavB and 1 kindL inside the cavA.
// Use CD-CBMC to build kindL

class IntraMoleculeExchange3 : public IntraMoleculeExchange1 {
public:
  IntraMoleculeExchange3(System &sys, StaticVals const &statV)
      : IntraMoleculeExchange1(sys, statV) {
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
  virtual double GetCoeff() const;
};

inline void IntraMoleculeExchange3::SetMEMC(StaticVals const &statV) {
  for (uint t = 0; t < exchangeRatioVec.size(); t++) {
    if (largeBBVec[t][0] != largeBBVec[t][1]) {
      printf("Error: In Intra-ME-3 move, two atoms with same name should be "
             "used as backbone.\n");
      printf("Atom names in backbone was set to %s or %s in %s residue.\n",
             statV.intraMemcVal.largeBBAtom1[t].c_str(),
             statV.intraMemcVal.largeBBAtom2[t].c_str(),
             statV.intraMemcVal.largeKind[t].c_str());
      exit(EXIT_FAILURE);
    }
  }
}

inline void IntraMoleculeExchange3::AdjustExRatio() {
  if (((counter + 1) % perAdjust) == 0) {
    int exMax = ceil((float)molInCavCount / (float)perAdjust);
    int exMin = 1;

    int index = kindS + kindL * molRef.GetKindsCount();
    double currAccept = (double)(accepted[sourceBox][index]) /
                        (double)(trial[sourceBox][index]);
    if (std::abs(currAccept - lastAccept) >= 0.05 * currAccept) {
      if (currAccept > lastAccept) {
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
    printf("ExchangeRatio: %d, Average kindS in cavity: %d \n", exchangeRatio,
           exMax);
  }
}

inline void IntraMoleculeExchange3::SetExchangeData() {
  IntraMoleculeExchange1::SetExchangeData();
}

inline uint IntraMoleculeExchange3::PickMolInCav() {
  uint state = mv::fail_state::NO_FAIL;
  // pick a random small kind in dense phase and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if (state == mv::fail_state::NO_FAIL) {
    centerA = comCurrRef.Get(pickedS);
    // Pick random vector and find two vectors that are perpendicular to V1
    XYZ cavVect = prng.RandomUnitVect();
    SetBasis(cavA, cavVect);
    // Calculate inverse matrix for cav here Inv = transpose
    TransposeMatrix(invCavA, cavA);

    // Find the small molecule kind in the cavityA
    if ((exchangeRatio == 1) ||
        calcEnRef.FindMolInCavity(molInCav, centerA, cavity, invCavA, sourceBox,
                                  kindS, exchangeRatio)) {
      // Find the exchangeRatio number of molecules kindS in cavity
      numInCavA = exchangeRatio;
      // add the random picked small molecule to the list.
      molIndexA.push_back(pickedS);
      kindIndexA.push_back(pickedKS);
      if (exchangeRatio == 1) {
        numSCavA = 1;
      } else {
        numSCavA = molInCav[kindS].size();
        // delete the picked small molecule from list
        for (uint s = 0; s < numSCavA; s++) {
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

      // pick a molecule from kindL in the same box
      numInCavB = 1;
      state = prng.PickMol(kindL, kindIndexB, molIndexB, numInCavB, sourceBox);
    } else {
      // reject the move
      state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    }

    // keep track of number of small molecule in cavity
    molInCavCount += numSCavA;
    counter++;
  }

  // After picking a large molecule, set the cavityB, and count kindS in cavityB
  if (state == mv::fail_state::NO_FAIL) {
    // use the predefine atom in kindL as the centerB
    uint start = molRef.MolStart(molIndexB[0]) + largeBB[0];
    centerB = coordCurrRef.Get(start);
    // Set V1 to a random vector and calculate two vector perpendicular to V1
    XYZ cavVect = prng.RandomUnitVect();
    SetBasis(cavB, cavVect);
    // Calculate inverse matrix for cav. Here Inv = Transpose
    TransposeMatrix(invCavB, cavB);
    if (exchangeRatio == 1) {
      numSCavB = 0;
    } else {
      // find how many of KindS exist in this centerB (COM of kindL)
      calcEnRef.FindMolInCavity(molInCav, centerB, cavity, invCavB, sourceBox,
                                kindS, exchangeRatio);
      numSCavB = molInCav[kindS].size();
    }
  }
  return state;
}

inline uint IntraMoleculeExchange3::Prep(const double subDraw,
                                         const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_INTRA_MEMC);
  // AdjustExRatio();
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    // transferring type A from source
    for (uint n = 0; n < numInCavA; n++) {
      newMolA.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef, sourceBox));
      oldMolA.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef, sourceBox));
    }

    for (uint n = 0; n < numInCavB; n++) {
      // transferring type B from source
      newMolB.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef, sourceBox));
      oldMolB.push_back(
          cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef, sourceBox));
    }

    // set the old coordinate after unwrap them
    for (uint n = 0; n < numInCavA; n++) {
      oldMolA[n].SetCoords(coordCurrRef, pStartA[n]);
      // copy cavA matrix to slant the old trial of molA
      oldMolA[n].SetCavMatrix(cavA);
      // set coordinate of moleA to newMolA, later it will shift to centerB
      newMolA[n].SetCoords(coordCurrRef, pStartA[n]);
      // copy cavB matrix to slant the new trial of molA
      newMolA[n].SetCavMatrix(cavB);
    }

    for (uint n = 0; n < numInCavB; n++) {
      oldMolB[n].SetCoords(coordCurrRef, pStartB[n]);
      // copy cavB matrix to slant the old trial of molB
      oldMolB[n].SetCavMatrix(cavB);
      // set coordinate of moleB to newMolB, later it will shift to centerA
      newMolB[n].SetCoords(coordCurrRef, pStartB[n]);
      // copy cavA matrix to slant the new trial of molB
      newMolB[n].SetCavMatrix(cavA);
    }

    // SetSeed(has cavity, COM is fixed, rotate around Backbone)
    for (uint n = 0; n < numInCavB; n++) {
      // Inserting molB from centerB to the centerA
      newMolB[n].SetSeed(centerA, cavity, true, true, true);
      // Set the Backbone of large molecule to be inserted
      newMolB[n].SetBackBone(largeBB);
      // perform rotational trial move for oldMolB
      oldMolB[n].SetSeed(centerB, cavity, true, true, true);
      // Set the Backbone of large molecule to be deleted
      oldMolB[n].SetBackBone(largeBB);
    }

    for (uint n = 0; n < numInCavA; n++) {
      if (n == 0) {
        // Inserting molA from cavity(centerA) to the cavityB(centerB)
        // COM is fixed, rotation around sphere
        newMolA[n].SetSeed(centerB, cavity, true, true, false);
        // perform trial move in cavity in sourceBox for oldMolA
        // COM is fixed, rotation around sphere
        oldMolA[n].SetSeed(centerA, cavity, true, true, false);
      } else {
        // Inserting molA from cavity(centerA) to the cavityB(centerB)
        newMolA[n].SetSeed(centerB, cavity, true, false, false);
        // perform trial move in cavity in sourceBox for oldMolA
        oldMolA[n].SetSeed(centerA, cavity, true, false, false);
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_INTRA_MEMC);
  return state;
}

inline uint IntraMoleculeExchange3::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_INTRA_MEMC);
  /// Remove the fixed COM kindS at the end because we insert it at first
  for (uint n = numInCavA; n > 0; --n) {
    cellList.RemoveMol(molIndexA[n - 1], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexA[n - 1]].BuildIDOld(oldMolA[n - 1],
                                               molIndexA[n - 1]);
    // Add bonded energy because we don't consider it in DCRotate.cpp
    oldMolA[n - 1].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n - 1]));
  }

  // Calc old energy before deleting
  for (uint n = 0; n < numInCavB; ++n) {
    cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexB[n]].BuildGrowOld(oldMolB[n], molIndexB[n]);
  }

  // Insert kindL to cavity of  center A using CD-CBMC
  for (uint n = 0; n < numInCavB; ++n) {
    molRef.kinds[kindIndexB[n]].BuildGrowNew(newMolB[n], molIndexB[n]);
    ShiftMol(n, false);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
    overlap |= newMolB[n].HasOverlap();
  }

  // Insert kindS to cavity of center B
  for (uint n = 0; n < numInCavA; ++n) {
    molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
    ShiftMol(n, true);
    cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
    // Add bonded energy because we don't consider it in DCRotate.cpp
    // If we find overlap, we still need to move the molecules so we can
    // reset things properly later, but we don't update the energy because
    // we will reject the move
    if (!overlap) {
      newMolA[n].AddEnergy(calcEnRef.MoleculeIntra(newMolA[n]));
      overlap |= newMolA[n].HasOverlap();
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_INTRA_MEMC);
  return mv::fail_state::NO_FAIL;
}

inline void IntraMoleculeExchange3::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_INTRA_MEMC);
  W_recip = 1.0;
  recipDiff = 0.0;
  correctDiff = 0.0;
  // No need to calculate the correction term for kindS since it is
  // inserted rigid body. We just need it for kindL
  if (!overlap) {
    for (uint n = 0; n < numInCavB; n++) {
      correctDiff += calcEwald->SwapCorrection(newMolB[n], molIndexB[n]);
      correctDiff -= calcEwald->SwapCorrection(oldMolB[n], molIndexB[n]);
    }
    // MolExchangeReciprocal returns the total change in recip energy. It
    // accumulates with each call, so we should use only the last of the two.
    recipDiff = calcEwald->MolExchangeReciprocal(newMolA, oldMolA, molIndexA,
                                                 molIndexA, true);
    recipDiff = calcEwald->MolExchangeReciprocal(newMolB, oldMolB, molIndexB,
                                                 molIndexB, false);

    W_recip = exp(-ffRef.beta * (recipDiff + correctDiff));
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_INTRA_MEMC);
}

inline double IntraMoleculeExchange3::GetCoeff() const {
  double ratioF = num::Factorial(numSCavA - 1) * num::Factorial(numSCavB) /
                  (num::Factorial(numSCavA - exchangeRatio) *
                   num::Factorial(numSCavB + exchangeRatio - 1));

  return ratioF;
}

inline void IntraMoleculeExchange3::Accept(const uint rejectState,
                                           const ulong step) {
  IntraMoleculeExchange1::Accept(rejectState, step);
}

#endif
