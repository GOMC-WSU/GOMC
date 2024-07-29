/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLCULETRANSFER_H
#define MOLCULETRANSFER_H

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC

#include "MoveBase.h"
#include "TrialMol.h"

//#define DEBUG_MOVES

class MoleculeTransfer : public MoveBase {
public:
  MoleculeTransfer(System &sys, StaticVals const &statV)
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
  double GetCoeff() const;
  uint GetBoxPairAndMol(const double subDraw, const double movPerc);
  MolPick molPick;
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  MoleculeLookup &molLookRef;
  Forcefield const &ffRef;
};

void MoleculeTransfer::PrintAcceptKind() {
  for (uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Mol-Transfer ",
           molRef.kinds[k].name.c_str());
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (moveSetRef.GetTrial(b, mv::MOL_TRANSFER, k) > 0)
        printf("%10.5f ",
               (100.0 * moveSetRef.GetAccept(b, mv::MOL_TRANSFER, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint MoleculeTransfer::GetBoxPairAndMol(const double subDraw,
                                               const double movPerc) {
  // Need to call a function to pick a molecule that is not fixed but cannot be
  // swap between boxes. (beta != 1, beta !=2)
  uint state = prng.PickMolAndBoxPair2(molIndex, kindIndex, sourceBox, destBox,
                                       subDraw, movPerc);
#if ENSEMBLE == GCMC
  if (state == mv::fail_state::NO_MOL_OF_KIND_IN_BOX && sourceBox == mv::BOX1) {
    std::cout << "Error: There are no molecules of kind "
              << molRef.kinds[kindIndex].name << " left in reservoir.\n";
    exit(EXIT_FAILURE);
  }
#endif

  if (state == mv::fail_state::NO_FAIL) {
    pStart = pLen = 0;
    molRef.GetRangeStartLength(pStart, pLen, molIndex);
  }
  return state;
}

inline uint MoleculeTransfer::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_SWAP);
  overlap = false;
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStart);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_SWAP);
  return state;
}

inline uint MoleculeTransfer::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_SWAP);
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);
  overlap = newMol.HasOverlap();
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_SWAP);
  return mv::fail_state::NO_FAIL;
}

inline void MoleculeTransfer::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_SWAP);
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
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_SWAP);
}

inline double MoleculeTransfer::GetCoeff() const {
#if ENSEMBLE == GEMC
  return (double)(molLookRef.NumKindInBoxSwappable(kindIndex, sourceBox)) /
         (double)(molLookRef.NumKindInBoxSwappable(kindIndex, destBox) + 1) *
         boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) { // Delete case
    if (ffRef.isFugacity) {
      return (double)(molLookRef.NumKindInBoxSwappable(kindIndex, sourceBox)) *
             boxDimRef.volInv[sourceBox] /
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return (double)(molLookRef.NumKindInBoxSwappable(kindIndex, sourceBox)) *
             boxDimRef.volInv[sourceBox] *
             exp(-BETA * molRef.kinds[kindIndex].chemPot);
    }
  } else { // Insertion case
    if (ffRef.isFugacity) {
      return boxDimRef.volume[destBox] /
             (double)(molLookRef.NumKindInBoxSwappable(kindIndex, destBox) +
                      1) *
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return boxDimRef.volume[destBox] /
             (double)(molLookRef.NumKindInBoxSwappable(kindIndex, destBox) +
                      1) *
             exp(BETA * molRef.kinds[kindIndex].chemPot);
    }
  }
#endif
}

inline void MoleculeTransfer::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_SWAP);
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
      // Add Reciprocal energy
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
      if (newMol.GetWeight() > SMALL_WEIGHT && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }

  } else // we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::MOL_TRANSFER, result, destBox, kindIndex);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_SWAP);
}

#endif

#endif
