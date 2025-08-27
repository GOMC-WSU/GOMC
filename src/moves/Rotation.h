/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef ROTATION_H
#define ROTATION_H

#include "MoveBase.h"

class Rotate : public MoveBase, public MolTransformBase {
public:
  Rotate(System &sys, StaticVals const &statV) : MoveBase(sys, statV) {}

  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const ulong step);
  virtual void PrintAcceptKind();

private:
  Intermolecular inter_LJ, inter_Real, recip;
};

void Rotate::PrintAcceptKind() {
  for (uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Rotation ", molRef.kinds[k].name.c_str());
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (moveSetRef.GetTrial(b, mv::ROTATE, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::ROTATE, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint Rotate::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_ROTATE);
  uint state = GetBoxAndMol(prng, molRef, subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL && molRef.NumAtoms(mk) <= 1)
    state = mv::fail_state::ROTATE_ON_SINGLE_ATOM;

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_ROTATE);
  return state;
}

inline uint Rotate::PrepNEMTMC(const uint box, const uint midx,
                               const uint kidx) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_ROTATE);
  b = box;
  m = midx;
  mk = kidx;
  pStart = pLen = 0;
  molRef.GetRangeStartLength(pStart, pLen, m);
  newMolPos.Uninit();
  newMolPos.Init(pLen);
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_ROTATE);
  return mv::fail_state::NO_FAIL;
}

inline uint Rotate::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_ROTATE);
  coordCurrRef.RotateRand(newMolPos, pStart, pLen, m, b,
                          moveSetRef.Scale(b, mv::ROTATE, mk));

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_ROTATE);
  return mv::fail_state::NO_FAIL;
}

inline void Rotate::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_ROTATE);
  cellList.RemoveMol(m, b, coordCurrRef);
  molRemoved = true;
  overlap = false;

  // calculate LJ interaction and real term of electrostatic interaction
  overlap = calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
  if (!overlap) {
    // calculate reciprocate term of electrostatic interaction
    recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_ROTATE);
}

inline void Rotate::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_ROTATE);
  bool res = false;

  if (rejectState == mv::fail_state::NO_FAIL) {
    double pr = prng();
    res =
        pr < exp(-BETA * (inter_LJ.energy + inter_Real.energy + recip.energy));
  }
  bool result = res && !overlap;

  if (result) {
    // Set new energy.
    // setting energy and virial of LJ interaction
    sysPotRef.boxEnergy[b].inter += inter_LJ.energy;
    // setting energy and virial of coulomb interaction
    sysPotRef.boxEnergy[b].real += inter_Real.energy;
    // setting energy and virial of recip term
    sysPotRef.boxEnergy[b].recip += recip.energy;

    // Copy coords
    newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
    if (recip.energy != 0.0) {
      calcEwald->UpdateRecip(b);
    }

    sysPotRef.Total();
    // Update the velocity
    velocity.UpdateMolVelocity(m, b);
  }

  if (molRemoved) {
    // It means that Recip energy is calculated and move not accepted
    if (!result && !overlap) {
      calcEwald->RestoreMol(m);
    }

    cellList.AddMol(m, b, coordCurrRef);
    molRemoved = false;
  }

  moveSetRef.Update(mv::ROTATE, result, b, mk);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_ROTATE);
}

#endif /*ROTATION_H*/
