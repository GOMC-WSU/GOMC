/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef TRANSLATE_H
#define TRANSLATE_H

#include "MoveBase.h"
#include "Rotation.h"

class Rotate;

class Translate : public MoveBase, public MolTransformBase {
public:
  Translate(System &sys, StaticVals const &statV) : MoveBase(sys, statV) {}

  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0);
  uint ReplaceRot(Rotate const &other);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint rejectState, const ulong step);
  virtual void PrintAcceptKind();

private:
  Intermolecular inter_LJ, inter_Real, recip;
  XYZ newCOM;
};

void Translate::PrintAcceptKind() {
  for (uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Displacement ",
           molRef.kinds[k].name.c_str());
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (moveSetRef.GetTrial(b, mv::DISPLACE, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::DISPLACE, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint Translate::ReplaceRot(Rotate const &other) {
  ReplaceWith(other);
  return mv::fail_state::NO_FAIL;
}

inline uint Translate::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_DISPLACE);
  uint state = GetBoxAndMol(prng, molRef, subDraw, movPerc);
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_DISPLACE);
  return state;
}

inline uint Translate::PrepNEMTMC(const uint box, const uint midx,
                                  const uint kidx) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_DISPLACE);
  b = box;
  m = midx;
  mk = kidx;
  pStart = pLen = 0;
  molRef.GetRangeStartLength(pStart, pLen, m);
  newMolPos.Uninit();
  newMolPos.Init(pLen);
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_DISPLACE);
  return mv::fail_state::NO_FAIL;
}

inline uint Translate::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_DISPLACE);
  coordCurrRef.TranslateRand(newMolPos, newCOM, pStart, pLen, m, b,
                             moveSetRef.Scale(b, mv::DISPLACE, mk));

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_DISPLACE);
  return mv::fail_state::NO_FAIL;
}

inline void Translate::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_DISPLACE);
  cellList.RemoveMol(m, b, coordCurrRef);
  molRemoved = true;
  overlap = false;

  // calculate LJ interaction and real term of electrostatic interaction
  overlap = calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
  if (!overlap) {
    // calculate reciprocate term of electrostatic interaction
    recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_DISPLACE);
}

inline void Translate::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_DISPLACE);
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
    ;

    // Copy coords
    newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
    comCurrRef.Set(m, newCOM);
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

  moveSetRef.Update(mv::DISPLACE, result, b, mk);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_DISPLACE);
}

#endif /*TRANSLATE_H*/
