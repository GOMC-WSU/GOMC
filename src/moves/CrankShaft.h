/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CRANKSHAFT_H
#define CRANKSHAFT_H

#include "MoveBase.h"
#include "TrialMol.h"

//#define DEBUG_MOVES

class CrankShaft : public MoveBase {
public:
  CrankShaft(System &sys, StaticVals const &statV)
      : MoveBase(sys, statV), molLookRef(sys.molLookupRef),
        ffRef(statV.forcefield) {}

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
  uint GetBoxAndMol(const double subDraw, const double movPerc);
  MolPick molPick;
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  double W_recip;
  double correct_old, correct_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular recipDiff;
  MoleculeLookup &molLookRef;
  Forcefield const &ffRef;
};

void CrankShaft::PrintAcceptKind() {
  for (uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Crank-Shaft ",
           molRef.kinds[k].name.c_str());
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (moveSetRef.GetTrial(b, mv::CRANKSHAFT, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::CRANKSHAFT, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint CrankShaft::GetBoxAndMol(const double subDraw,
                                     const double movPerc) {
#if ENSEMBLE == GCMC
  sourceBox = mv::BOX0;
  uint state = prng.PickMol(molIndex, kindIndex, sourceBox, subDraw, movPerc);
#else
  uint state =
      prng.PickMolAndBox(molIndex, kindIndex, sourceBox, subDraw, movPerc);
#endif

  // molecule will be removed and insert in same box
  destBox = sourceBox;

  if (state == mv::fail_state::NO_FAIL) {
    pStart = pLen = 0;
    molRef.GetRangeStartLength(pStart, pLen, molIndex);
  }
  return state;
}

inline uint CrankShaft::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_CRANKSHAFT);
  overlap = false;
  uint state = GetBoxAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStart);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_CRANKSHAFT);
  return state;
}

inline uint CrankShaft::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_CRANKSHAFT);
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  molRef.kinds[kindIndex].CrankShaft(oldMol, newMol, molIndex);
  overlap = newMol.HasOverlap();
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_CRANKSHAFT);
  return mv::fail_state::NO_FAIL;
}

inline void CrankShaft::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_CRANKSHAFT);
  // since number of molecules would not change in the box,
  W_recip = 1.0;
  correct_old = 0.0;
  correct_new = 0.0;

  if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
    correct_new = calcEwald->SwapCorrection(newMol, molIndex);
    correct_old = calcEwald->SwapCorrection(oldMol, molIndex);
    recipDiff.energy =
        calcEwald->MolReciprocal(newMol.GetCoords(), molIndex, sourceBox);
    // self energy is same
    W_recip =
        exp(-1.0 * ffRef.beta * (recipDiff.energy + correct_new - correct_old));
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_CRANKSHAFT);
}

inline void CrankShaft::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_CRANKSHAFT);
  bool result;
  // If we didn't skip the move calculation
  if (rejectState == mv::fail_state::NO_FAIL) {
    double Wo = oldMol.GetWeight();
    double Wn = newMol.GetWeight();
    double Wrat = Wn / Wo * W_recip;

    // safety to make sure move will be rejected in overlap case
    if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
      result = prng() < Wrat;
    } else
      result = false;

    if (result) {
      // Add rest of energy.
      sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
      sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
      //Add Reciprocal energy difference
      sysPotRef.boxEnergy[destBox].recip += recipDiff.energy;
      // Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
      sysPotRef.boxEnergy[destBox].correction += correct_new;

      // Set coordinates, new COM; shift index to new box's list
      newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
      comCurrRef.SetNew(molIndex, destBox);
      cellList.AddMol(molIndex, destBox, coordCurrRef);

      // Zero out box energies to prevent small number
      // errors in double.
      if (molLookRef.NumInBox(sourceBox) == 1) {
        sysPotRef.boxEnergy[sourceBox].inter = 0;
        sysPotRef.boxVirial[sourceBox].inter = 0;
        sysPotRef.boxEnergy[sourceBox].real = 0;
        sysPotRef.boxVirial[sourceBox].real = 0;
      }

      calcEwald->UpdateRecip(sourceBox);

      // Retotal
      sysPotRef.Total();
      // Update the velocity
      velocity.UpdateMolVelocity(molIndex, sourceBox);
    } else {
      cellList.AddMol(molIndex, sourceBox, coordCurrRef);

      // when weight is 0, MolDestSwap() will not be executed, thus cos/sin
      // molRef will not be changed. Also since no memcpy, doing restore
      // results in memory overwrite
      if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }
  } else // else we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::CRANKSHAFT, result, sourceBox, kindIndex);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_CRANKSHAFT);
}

#endif
