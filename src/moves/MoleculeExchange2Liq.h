/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLECULEEXCHANGE2LIQ_H
#define MOLECULEEXCHANGE2LIQ_H

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC

#include "GeomLib.h"
#include "MoleculeExchange1.h"
#include "TrialMol.h"
#include <cmath>

using std::vector;
using namespace geom;

// MEMC-2-Liq Move:
//
// Swapping one Large molecule with one or more small molecules and vice versa.
// Sub-Volume location and orientation is based on the COM and backbone of the
// the small and large molecule. We always pick small molecule from sourceBox
// and Large from destBox to swap, and insert in the cavity left behind by small
// or large molecule kind The acceptance criteria for GCMC simulation is similar
// to MEMC-2, but the GEMC simulation has different acceptance criteria

class MoleculeExchange2Liq : public MoleculeExchange1 {
public:
  MoleculeExchange2Liq(System &sys, StaticVals const &statV)
      : MoleculeExchange1(sys, statV) {
    if (enableID) {
      SetMEMC(statV);
    }
    for (int b = 0; b < BOX_TOTAL; ++b) {
      cavityMatrix[b].Init(3);
      invCavityMatrix[b].Init(3);
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
  virtual void SetMEMC(StaticVals const &statV);
  virtual void SetExchangeData();
  virtual uint PickMolInCav();
  virtual double GetCoeff() const;

  uint SetBoxMolKind(const double subDraw, const double movPerc);
  void ShiftMol(uint box, const uint n, const uint from, const uint to);
  void RecoverMol(uint box, const uint n, const uint from, const uint to);
  int smallBB[2];

  vector<uint> pStart[BOX_TOTAL], pLen[BOX_TOTAL];
  vector<uint> molIndex[BOX_TOTAL], kindIndex[BOX_TOTAL];
  vector<cbmc::TrialMol> oldMol[BOX_TOTAL], newMol[BOX_TOTAL];
  // To keep track of small molecules in cavity in each box
  uint numMolInCavity[BOX_TOTAL], totalMolInCavity[BOX_TOTAL];
  // To store total sets of exchange pairs
  vector<vector<uint>> smallBBVec;
  vector<vector<uint>> moleculeIndexInCavity[BOX_TOTAL];
  XYZ cavityCenter[BOX_TOTAL];
  XYZArray cavityMatrix[BOX_TOTAL], invCavityMatrix[BOX_TOTAL];
};

inline void MoleculeExchange2Liq::SetMEMC(StaticVals const &statV) {
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
        printf("Error: In ME-2-Liq move, atom name %s or %s was not found in "
               "%s residue.\n",
               statV.memcVal.smallBBAtom1[t].c_str(),
               statV.memcVal.smallBBAtom2[t].c_str(),
               statV.memcVal.smallKind[t].c_str());
        exit(EXIT_FAILURE);
      }
    }

    if (molRef.kinds[kindSVec[t]].NumAtoms() > 1) {
      if (smallBB[0] == smallBB[1]) {
        printf("Error: In ME-2-Liq move, atom names in small molecule backbone "
               "cannot be same!\n");
        exit(EXIT_FAILURE);
      }
    }
    vector<uint> temp(smallBB, smallBB + 2);
    smallBBVec.push_back(temp);
  }
}

inline uint MoleculeExchange2Liq::SetBoxMolKind(const double subDraw,
                                                const double movPerc) {
  uint state = mv::fail_state::NO_FAIL;
  overlap = false;
  // pick a random source and dest Box. We always replace Large kind in
  //  dest box with small kind in source box
  prng.PickBoxPair(sourceBox, destBox, subDraw, movPerc);
  // pick one of the exchange type
  SetExchangeData();

  // adjust exchange rate based on number of small kind in cavity
  // AdjustExRatio();

  for (uint b = 0; b < BOX_TOTAL; ++b) {
    molIndex[b].clear();
    kindIndex[b].clear();
    oldMol[b].clear();
    newMol[b].clear();
  }

  // Pick molecules of kind S and L in each box, setup the cavity matrix
  state = PickMolInCav();

  if (state == mv::fail_state::NO_FAIL) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      pStart[b].clear();
      pLen[b].clear();
      pStart[b].resize(numMolInCavity[b]);
      pLen[b].resize(numMolInCavity[b]);

      for (uint n = 0; n < numMolInCavity[b]; n++) {
        pStart[b][n] = pLen[b][n] = 0;
        molRef.GetRangeStartLength(pStart[b][n], pLen[b][n], molIndex[b][n]);
      }
    }
  }

  return state;
}

inline void MoleculeExchange2Liq::SetExchangeData() {
  uint exType = prng.randIntExc(exchangeRatioVec.size());
  kindS = kindSVec[exType];
  kindL = kindLVec[exType];
  exchangeRatio = exchangeRatioVec[exType];
  largeBB[0] = largeBBVec[exType][0];
  largeBB[1] = largeBBVec[exType][1];
  smallBB[0] = smallBBVec[exType][0];
  smallBB[1] = smallBBVec[exType][1];
}

inline uint MoleculeExchange2Liq::PickMolInCav() {
  uint state = mv::fail_state::NO_FAIL;
  // Exchange N_ex small molecule with one large kind
  numMolInCavity[sourceBox] = exchangeRatio;
  numMolInCavity[destBox] = 1;
  // pick a random small kind in source box and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if (state == mv::fail_state::NO_FAIL) {
    cavityCenter[sourceBox] = comCurrRef.Get(pickedS);
    if (molRef.NumAtoms(kindS) == 1) {
      // with single site molecules, we cannot define the backbone, therefore
      // use random orientation
      SetBasis(cavityMatrix[sourceBox], prng.RandomUnitVect());
    } else {
      // Orient the cavity with backbone of picked small molecule
      uint start = molRef.MolStart(pickedS) + smallBB[0];
      uint end = molRef.MolStart(pickedS) + smallBB[1];
      SetBasis(
          cavityMatrix[sourceBox],
          boxDimRef.MinImage(coordCurrRef.Difference(start, end), sourceBox));
    }
    // Calculate inverse matrix for cav here Inv = transpose
    TransposeMatrix(invCavityMatrix[sourceBox], cavityMatrix[sourceBox]);

    // Find all the molecule small kind in the cavity of sourceBox. Returns
    // False if we cannot find exchangeRatio of small molecule in source box
    if ((exchangeRatio == 1) ||
        calcEnRef.FindMolInCavity(
            moleculeIndexInCavity[sourceBox], cavityCenter[sourceBox], cavity,
            invCavityMatrix[sourceBox], sourceBox, kindS, exchangeRatio)) {
      // Find the exchangeRatio number of small molecules kind in cavity
      // add the random picked small molecule to the list.
      molIndex[sourceBox].push_back(pickedS);
      kindIndex[sourceBox].push_back(kindS);
      if (exchangeRatio == 1) {
        totalMolInCavity[sourceBox] = 1;
      } else {
        totalMolInCavity[sourceBox] =
            moleculeIndexInCavity[sourceBox][kindS].size();
        // delete the picked small molecule from list
        for (uint s = 0; s < totalMolInCavity[sourceBox]; s++) {
          if (pickedS == moleculeIndexInCavity[sourceBox][kindS][s])
            moleculeIndexInCavity[sourceBox][kindS].erase(
                moleculeIndexInCavity[sourceBox][kindS].begin() + s);
        }
      }

      for (uint n = 1; n < numMolInCavity[sourceBox]; n++) {
        // pick random exchangeRatio number of kindS in cavity
        uint picked =
            prng.randIntExc(moleculeIndexInCavity[sourceBox][kindS].size());
        molIndex[sourceBox].push_back(
            moleculeIndexInCavity[sourceBox][kindS][picked]);
        kindIndex[sourceBox].push_back(
            molRef.GetMolKind(molIndex[sourceBox][n]));
        moleculeIndexInCavity[sourceBox][kindS].erase(
            moleculeIndexInCavity[sourceBox][kindS].begin() + picked);
      }

      // pick a Large molecule kind from destBox
      uint pickedL, pickedKL;
      state = prng.PickMol(kindL, pickedKL, pickedL, destBox);
      // Now set up the cavity for large molecule in destBox
      if (state == mv::fail_state::NO_FAIL) {
        molIndex[destBox].push_back(pickedL);
        kindIndex[destBox].push_back(kindL);
        // Find the center of large molecule in destbox
        cavityCenter[destBox] = comCurrRef.Get(pickedL);
        // If we want to orient the cavity with backbone of picked large mol
        if (molRef.NumAtoms(kindL) == 1) {
          // with single site molecules, we cannot define the backbone,
          // therefore use random orientation
          SetBasis(cavityMatrix[destBox], prng.RandomUnitVect());
        } else {
          // Orient the cavity with backbone of picked large molecule kind
          uint start = molRef.MolStart(pickedL) + largeBB[0];
          uint end = molRef.MolStart(pickedL) + largeBB[1];
          SetBasis(
              cavityMatrix[destBox],
              boxDimRef.MinImage(coordCurrRef.Difference(start, end), destBox));
        }
        // Calculate inverse matrix for cav here Inv = transpose
        TransposeMatrix(invCavityMatrix[destBox], cavityMatrix[destBox]);

        if (exchangeRatio == 1) {
          totalMolInCavity[destBox] = 0;
        } else {
          // find how many of KindS exist in the cavity of destBox
          calcEnRef.FindMolInCavity(
              moleculeIndexInCavity[destBox], cavityCenter[destBox], cavity,
              invCavityMatrix[destBox], destBox, kindS, exchangeRatio);
          totalMolInCavity[destBox] =
              moleculeIndexInCavity[destBox][kindS].size();
        }
      }
    } else {
      // reject the move
      state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    }
  }

  return state;
}

inline uint MoleculeExchange2Liq::Prep(const double subDraw,
                                       const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MEMC);
  uint state = SetBoxMolKind(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    numTypeASource = (double)(molLookRef.NumKindInBoxSwappable(
        kindIndex[sourceBox][0], sourceBox)); // small kind in sourceBox
    numTypeADest = (double)(molLookRef.NumKindInBoxSwappable(
        kindIndex[sourceBox][0], destBox)); // small kind in destBox
    numTypeBSource = (double)(molLookRef.NumKindInBoxSwappable(
        kindIndex[destBox][0], sourceBox)); // large kind in source box
    numTypeBDest = (double)(molLookRef.NumKindInBoxSwappable(
        kindIndex[destBox][0], destBox)); // large kind in destBox

    // transferring small kind from sourceBox to destBox
    for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
      newMol[sourceBox].push_back(cbmc::TrialMol(
          molRef.kinds[kindIndex[sourceBox][n]], boxDimRef, destBox));
      oldMol[sourceBox].push_back(cbmc::TrialMol(
          molRef.kinds[kindIndex[sourceBox][n]], boxDimRef, sourceBox));
    }

    // transferring large kind from destBox to sourceBox
    for (uint n = 0; n < numMolInCavity[destBox]; n++) {
      newMol[destBox].push_back(cbmc::TrialMol(
          molRef.kinds[kindIndex[destBox][n]], boxDimRef, sourceBox));
      oldMol[destBox].push_back(cbmc::TrialMol(
          molRef.kinds[kindIndex[destBox][n]], boxDimRef, destBox));
    }

    // set the oldMol and newMol coordinate for small kind, after proper wrap &
    // unwrap
    for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
      XYZArray mol(pLen[sourceBox][n]);
      coordCurrRef.CopyRange(mol, pStart[sourceBox][n], 0, pLen[sourceBox][n]);
      boxDimRef.UnwrapPBC(mol, sourceBox,
                          comCurrRef.Get(molIndex[sourceBox][n]));
      boxDimRef.WrapPBC(mol, destBox);
      oldMol[sourceBox][n].SetCoords(coordCurrRef, pStart[sourceBox][n]);
      // set coordinate of oldMol to newMol for small kind, later it will shift
      // to cavityCenter
      newMol[sourceBox][n].SetCoords(mol, 0);
      // copy cavityMatrix of small kind in sourceBox for the old trial
      oldMol[sourceBox][n].SetCavMatrix(cavityMatrix[sourceBox]);
      // copy cavityMatrix of large kind in dest for the new trial
      newMol[sourceBox][n].SetCavMatrix(cavityMatrix[destBox]);
    }

    // set the oldMol and newMol coordinate for large kind, after proper wrap &
    // unwrap
    for (uint n = 0; n < numMolInCavity[destBox]; n++) {
      XYZArray mol(pLen[destBox][n]);
      coordCurrRef.CopyRange(mol, pStart[destBox][n], 0, pLen[destBox][n]);
      boxDimRef.UnwrapPBC(mol, destBox, comCurrRef.Get(molIndex[destBox][n]));
      boxDimRef.WrapPBC(mol, sourceBox);
      oldMol[destBox][n].SetCoords(coordCurrRef, pStart[destBox][n]);
      // set coordinate of oldMol to newMol for large kind, later it will shift
      // to cavityCenter
      newMol[destBox][n].SetCoords(mol, 0);
      // copy cavityMatrix of large kind in destBox for the old trial
      oldMol[destBox][n].SetCavMatrix(cavityMatrix[destBox]);
      // copy cavityMatrix of small kind in source for the new trial
      newMol[destBox][n].SetCavMatrix(cavityMatrix[sourceBox]);
    }

    // set cavity center and dimension for large kind
    for (uint n = 0; n < numMolInCavity[destBox]; n++) {
      // SetSeed(x, x, has cavity, COM is fixed, rotate around Backbone)

      // Inserting Large kind from destBox to the cavity in sourceBox
      newMol[destBox][n].SetSeed(cavityCenter[sourceBox], cavity, true, true,
                                 true);
      // Set the Backbone of large molecule to be inserted
      newMol[destBox][n].SetBackBone(largeBB);
      // perform rotational trial move around backbone for large kind in destBox
      oldMol[destBox][n].SetSeed(cavityCenter[destBox], cavity, true, true,
                                 true);
      // Set the Backbone of large molecule to perform rotational trial around
      // it
      oldMol[destBox][n].SetBackBone(largeBB);
    }

    // set cavity center and dimension for small large kind
    for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
      // SetSeed(x, x, has cavity, COM is fixed, rotate around Backbone)

      if (n == 0) {
        // Inserting Small kind from source to the cavity in destBox
        newMol[sourceBox][n].SetSeed(cavityCenter[destBox], cavity, true, true,
                                     true);
        // Set the Backbone of small molecule to be inserted
        newMol[sourceBox][n].SetBackBone(smallBB);
        // Perform rotational trial move around backbone for small kind in
        // sourceBox
        oldMol[sourceBox][n].SetSeed(cavityCenter[sourceBox], cavity, true,
                                     true, true);
        // Set the Backbone of small molecule to perform rotational trial around
        // it
        oldMol[sourceBox][n].SetBackBone(smallBB);
      } else {
        // Inserting Small kind from sourceBox to the cavity in destBox
        newMol[sourceBox][n].SetSeed(cavityCenter[destBox], cavity, true, false,
                                     false);
        // Perform trial move for small kind in cavity in sourceBox
        oldMol[sourceBox][n].SetSeed(cavityCenter[sourceBox], cavity, true,
                                     false, false);
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MEMC);
  return state;
}

inline uint MoleculeExchange2Liq::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MEMC);

  // Calculate tail correction
  W_tc = 1.0;
  if (ffRef.useLRC) {
    double delTC = 0.0;
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      uint *kCount = new uint[molRef.kindsCount];
      for (uint k = 0; k < molRef.kindsCount; ++k) {
        kCount[k] = molLookRef.NumKindInBox(k, b);
      }

      // Moving small kind from sourceBox to destBox
      // Moving large kind from destBox to sourceBox
      if (b == sourceBox) {
        kCount[kindIndex[sourceBox][0]] -= numMolInCavity[sourceBox];
        kCount[kindIndex[destBox][0]] += numMolInCavity[destBox];
      } else if (b == destBox) {
        kCount[kindIndex[sourceBox][0]] += numMolInCavity[sourceBox];
        kCount[kindIndex[destBox][0]] -= numMolInCavity[destBox];
      }
      tcNew[b].energy = calcEnRef.EnergyCorrection(b, kCount);
      delTC += tcNew[b].energy - sysPotRef.boxEnergy[b].tailCorrection;
      delete[] kCount;
    }
    W_tc = std::exp(-ffRef.beta * delTC);
  }

  // Calc Old energy of small molecule and delete it from sourceBox
  // Remove the fixed COM small mol at the end because we insert it at first
  for (uint n = numMolInCavity[sourceBox]; n > 0; --n) {
    cellList.RemoveMol(molIndex[sourceBox][n - 1], sourceBox, coordCurrRef);
    molRef.kinds[kindIndex[sourceBox][n - 1]].BuildIDOld(
        oldMol[sourceBox][n - 1], molIndex[sourceBox][n - 1]);
    // If we find overlap, we still need to move the molecules so we can
    // reset things properly later, but we don't update the energy because
    // we will reject the move
    overlap |= oldMol[sourceBox][n - 1].HasOverlap();
    if (!overlap) {
      // Add bonded energy because we don't consider it in DCRotate.cpp
      oldMol[sourceBox][n - 1].AddEnergy(
          calcEnRef.MoleculeIntra(oldMol[sourceBox][n - 1]));
    }
  }

  // Calc Old energy of large molecule and delete it from destBox
  for (uint n = 0; n < numMolInCavity[destBox]; ++n) {
    cellList.RemoveMol(molIndex[destBox][n], destBox, coordCurrRef);
    molRef.kinds[kindIndex[destBox][n]].BuildIDOld(oldMol[destBox][n],
                                                   molIndex[destBox][n]);
    // If we find overlap, we still need to move the molecules so we can
    // reset things properly later, but we don't update the energy because
    // we will reject the move
    overlap |= oldMol[destBox][n].HasOverlap();
    if (!overlap) {
      oldMol[destBox][n].AddEnergy(calcEnRef.MoleculeIntra(oldMol[destBox][n]));
    }
  }

  // Insert small molecule to destBox and calc new energy
  for (uint n = 0; n < numMolInCavity[sourceBox]; ++n) {
    molRef.kinds[kindIndex[sourceBox][n]].BuildIDNew(newMol[sourceBox][n],
                                                     molIndex[sourceBox][n]);
    this->ShiftMol(sourceBox, n, sourceBox, destBox);
    cellList.AddMol(molIndex[sourceBox][n], destBox, coordCurrRef);
    // If we find overlap, we still need to move the molecules so we can
    // reset things properly later, but we don't update the energy because
    // we will reject the move
    overlap |= newMol[sourceBox][n].HasOverlap();
    if (!overlap) {
      // Add bonded energy because we don't consider it in DCRotate.cpp
      newMol[sourceBox][n].AddEnergy(
          calcEnRef.MoleculeIntra(newMol[sourceBox][n]));
    }
  }

  // Insert large molecule to sourceBox and calc new energy
  for (uint n = 0; n < numMolInCavity[destBox]; ++n) {
    molRef.kinds[kindIndex[destBox][n]].BuildIDNew(newMol[destBox][n],
                                                   molIndex[destBox][n]);
    this->ShiftMol(destBox, n, destBox, sourceBox);
    cellList.AddMol(molIndex[destBox][n], sourceBox, coordCurrRef);
    // If we find overlap, we still need to move the molecules so we can
    // reset things properly later, but we don't update the energy because
    // we will reject the move
    overlap |= newMol[destBox][n].HasOverlap();
    if (!overlap) {
      // Add bonded energy because we don't consider it in DCRotate.cpp
      newMol[destBox][n].AddEnergy(calcEnRef.MoleculeIntra(newMol[destBox][n]));
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MEMC);
  return mv::fail_state::NO_FAIL;
}

inline void MoleculeExchange2Liq::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MEMC);
  W_recip = 1.0;
  correct_oldA = 0.0, correct_newA = 0.0;
  self_oldA = 0.0, self_newA = 0.0;
  correct_oldB = 0.0, correct_newB = 0.0;
  self_oldB = 0.0, self_newB = 0.0;
  recipDest = 0.0, recipSource = 0.0;

  if (!overlap) {
    for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
      correct_newA += calcEwald->SwapCorrection(newMol[sourceBox][n]);
      correct_oldA += calcEwald->SwapCorrection(oldMol[sourceBox][n]);
      self_newA += calcEwald->SwapSelf(newMol[sourceBox][n]);
      self_oldA += calcEwald->SwapSelf(oldMol[sourceBox][n]);
    }
    recipDest = calcEwald->MolExchangeReciprocal(
        newMol[sourceBox], oldMol[destBox], molIndex[sourceBox],
        molIndex[destBox], true);

    for (uint n = 0; n < numMolInCavity[destBox]; n++) {
      correct_newB += calcEwald->SwapCorrection(newMol[destBox][n]);
      correct_oldB += calcEwald->SwapCorrection(oldMol[destBox][n]);
      self_newB += calcEwald->SwapSelf(newMol[destBox][n]);
      self_oldB += calcEwald->SwapSelf(oldMol[destBox][n]);
    }
    recipSource = calcEwald->MolExchangeReciprocal(
        newMol[destBox], oldMol[sourceBox], molIndex[destBox],
        molIndex[sourceBox], true);

    // need to contribute the self and correction energy
    W_recip =
        std::exp(-ffRef.beta * (recipSource + recipDest + correct_newA -
                                correct_oldA + correct_newB - correct_oldB +
                                self_newA - self_oldA + self_newB - self_oldB));
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MEMC);
}

inline double MoleculeExchange2Liq::GetCoeff() const {
#if ENSEMBLE == GEMC
  double ratioM = numTypeASource * numTypeBDest /
                  ((numTypeBSource + 1.0) * (numTypeADest + exchangeRatio));
  double ratioF =
      num::Factorial(totalMolInCavity[sourceBox] - 1) /
      (num::Factorial(totalMolInCavity[destBox], exchangeRatio - 1) *
       num::Factorial(totalMolInCavity[sourceBox] - exchangeRatio));

  // subVolume will cancel out
  return ratioF * ratioM;
#elif ENSEMBLE == GCMC
  if (ffRef.isFugacity) {
    if (sourceBox == mv::BOX0) {
      // Insert Large molecule
      double delA = molRef.kinds[kindIndex[sourceBox][0]].chemPot *
                    numMolInCavity[sourceBox];
      double insB =
          molRef.kinds[kindIndex[destBox][0]].chemPot * numMolInCavity[destBox];
      double ratioF =
          num::Factorial(totalMolInCavity[sourceBox] - 1) /
          num::Factorial(totalMolInCavity[sourceBox] - exchangeRatio);
      double ratioM = numTypeASource / (numTypeBSource + 1.0);
      return (insB / delA) * ratioF * ratioM / pow(volCav, exchangeRatio - 1);
    } else {
      // Delete Large Molecule
      double insA = molRef.kinds[kindIndex[sourceBox][0]].chemPot *
                    numMolInCavity[sourceBox];
      double delB =
          molRef.kinds[kindIndex[destBox][0]].chemPot * numMolInCavity[destBox];
      double ratioF =
          num::Factorial(totalMolInCavity[destBox]) /
          num::Factorial(totalMolInCavity[destBox] + exchangeRatio - 1);
      double ratioM = numTypeBDest / (numTypeADest + exchangeRatio);
      return (insA / delB) * ratioF * ratioM * pow(volCav, exchangeRatio - 1);
    }
  } else {
    if (sourceBox == mv::BOX0) {
      // Insert Large molecule
      double delA = (-BETA * molRef.kinds[kindIndex[sourceBox][0]].chemPot *
                     numMolInCavity[sourceBox]);
      double insB = (BETA * molRef.kinds[kindIndex[destBox][0]].chemPot *
                     numMolInCavity[destBox]);
      double ratioF =
          num::Factorial(totalMolInCavity[sourceBox] - 1) /
          num::Factorial(totalMolInCavity[sourceBox] - exchangeRatio);
      double ratioM = numTypeASource / (numTypeBSource + 1.0);
      return std::exp(delA + insB) * ratioF * ratioM /
             pow(volCav, exchangeRatio - 1);
    } else {
      // Delete Large Molecule
      double insA = (BETA * molRef.kinds[kindIndex[sourceBox][0]].chemPot *
                     numMolInCavity[sourceBox]);
      double delB = (-BETA * molRef.kinds[kindIndex[destBox][0]].chemPot *
                     numMolInCavity[destBox]);
      double ratioF =
          num::Factorial(totalMolInCavity[destBox]) /
          num::Factorial(totalMolInCavity[destBox] + exchangeRatio - 1);
      double ratioM = numTypeBDest / (numTypeADest + exchangeRatio);
      return std::exp(insA + delB) * ratioF * ratioM *
             pow(volCav, exchangeRatio - 1);
    }
  }
#endif
}

inline void MoleculeExchange2Liq::ShiftMol(uint box, const uint n,
                                           const uint from, const uint to) {

  // Shift kind in box, from to box
  newMol[box][n].GetCoords().CopyRange(coordCurrRef, 0, pStart[box][n],
                                       pLen[box][n]);
  comCurrRef.SetNew(molIndex[box][n], to);
  molLookRef.ShiftMolBox(molIndex[box][n], from, to, kindIndex[box][n]);
}

inline void MoleculeExchange2Liq::RecoverMol(uint box, const uint n,
                                             const uint from, const uint to) {
  XYZArray mol(pLen[box][n]);
  oldMol[box][n].GetCoords().CopyRange(mol, 0, 0, pLen[box][n]);
  boxDimRef.WrapPBC(mol, to);
  mol.CopyRange(coordCurrRef, 0, pStart[box][n], pLen[box][n]);
  comCurrRef.SetNew(molIndex[box][n], to);
  molLookRef.ShiftMolBox(molIndex[box][n], from, to, kindIndex[box][n]);
}

inline void MoleculeExchange2Liq::Accept(const uint rejectState,
                                         const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_MEMC);
  bool result = !overlap;

  // If we didn't skip the move calculation
  if (rejectState == mv::fail_state::NO_FAIL) {
    // If there is no overlap then calculate if we should accept the move based
    // on the coefficient and our random value
    if (result) {
      double molTransCoeff = GetCoeff();
      double Wrat = W_tc * W_recip;

      for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
        Wrat *=
            newMol[sourceBox][n].GetWeight() / oldMol[sourceBox][n].GetWeight();
      }

      for (uint n = 0; n < numMolInCavity[destBox]; n++) {
        Wrat *= newMol[destBox][n].GetWeight() / oldMol[destBox][n].GetWeight();
      }

      result = prng() < molTransCoeff * Wrat;
    }

    if (result) {
      // Add tail corrections
      sysPotRef.boxEnergy[sourceBox].tailCorrection = tcNew[sourceBox].energy;
      sysPotRef.boxEnergy[destBox].tailCorrection = tcNew[destBox].energy;

      // Add rest of energy.
      for (uint n = 0; n < numMolInCavity[destBox]; n++) {
        sysPotRef.boxEnergy[sourceBox] += newMol[destBox][n].GetEnergy();
        sysPotRef.boxEnergy[destBox] -= oldMol[destBox][n].GetEnergy();
      }

      for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
        sysPotRef.boxEnergy[sourceBox] -= oldMol[sourceBox][n].GetEnergy();
        sysPotRef.boxEnergy[destBox] += newMol[sourceBox][n].GetEnergy();
      }

      // Add Reciprocal energy
      sysPotRef.boxEnergy[sourceBox].recip += recipSource;
      sysPotRef.boxEnergy[destBox].recip += recipDest;
      // Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_oldA;
      sysPotRef.boxEnergy[sourceBox].correction += correct_newB;
      sysPotRef.boxEnergy[destBox].correction += correct_newA;
      sysPotRef.boxEnergy[destBox].correction -= correct_oldB;
      // Add self energy
      sysPotRef.boxEnergy[sourceBox].self -= self_oldA;
      sysPotRef.boxEnergy[sourceBox].self += self_newB;
      sysPotRef.boxEnergy[destBox].self += self_newA;
      sysPotRef.boxEnergy[destBox].self -= self_oldB;

      // If recip energy is unchanged, the SumI and SumR arrays are unchanged
      if (recipSource != 0.0 || recipDest != 0.0) {
        calcEwald->UpdateRecip(sourceBox);
        calcEwald->UpdateRecip(destBox);
      }

      // small and large molecule have already been transferred to destBox and
      // added to cellList, so don't need to update the cellList

      // Recalculate total
      sysPotRef.Total();

      // Update the velocity
      for (uint n = 0; n < numMolInCavity[destBox]; n++) {
        velocity.UpdateMolVelocity(molIndex[destBox][n], sourceBox);
      }
      for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
        velocity.UpdateMolVelocity(molIndex[sourceBox][n], destBox);
      }

    } else {
      // transfer back small molecule kind from destBox to source
      for (uint n = 0; n < numMolInCavity[sourceBox]; n++) {
        cellList.RemoveMol(molIndex[sourceBox][n], destBox, coordCurrRef);
        this->RecoverMol(sourceBox, n, destBox, sourceBox);
        cellList.AddMol(molIndex[sourceBox][n], sourceBox, coordCurrRef);
      }
      // transfer back large molecule kind from sourceBox to dest
      for (uint n = 0; n < numMolInCavity[destBox]; n++) {
        cellList.RemoveMol(molIndex[destBox][n], sourceBox, coordCurrRef);
        this->RecoverMol(destBox, n, sourceBox, destBox);
        cellList.AddMol(molIndex[destBox][n], destBox, coordCurrRef);
      }
    }
  } else {
    // else we didn't even try because we knew it would fail
    result = false;
  }

  moveSetRef.Update(mv::MEMC, result, sourceBox);
  moveSetRef.Update(mv::MEMC, result, destBox);

  // If we consider total acceptance of S->L and L->S
  AcceptKind(result, kindS + kindL * molRef.GetKindsCount(), sourceBox);
  AcceptKind(result, kindS + kindL * molRef.GetKindsCount(), destBox);

  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_MEMC);
}

#endif

#endif
