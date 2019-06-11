/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ROTATION_H
#define ROTATION_H

#include "MoveBase.h"


class Rotate : public MoveBase, public MolTransformBase
{
public:
  Rotate(System &sys, StaticVals const& statV) : MoveBase(sys, statV) {}

  virtual uint Prep(const real subDraw, const real movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);
  virtual void PrintAcceptKind();
private:
  Intermolecular inter_LJ, inter_Real, recip;
};

void Rotate::PrintAcceptKind()
{
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Rotation ", molRef.kinds[k].name.c_str());
    for(uint b = 0; b < BOX_TOTAL; b++) {
      if(moveSetRef.GetTrial(b, mv::ROTATE, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::ROTATE, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint Rotate::Prep(const real subDraw, const real movPerc)
{
  uint state = GetBoxAndMol(prng, molRef, subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL && molRef.NumAtoms(mk)  <= 1)
    state = mv::fail_state::ROTATE_ON_SINGLE_ATOM;
  return state;
}

inline uint Rotate::Transform()
{
  coordCurrRef.RotateRand(newMolPos, pStart, pLen, m, b,
                          moveSetRef.Scale(b, mv::ROTATE, mk));
  return mv::fail_state::NO_FAIL;
}

inline void Rotate::CalcEn()
{
  cellList.RemoveMol(m, b, coordCurrRef);
  molRemoved = true;
  overlap = false;

  //calculate LJ interaction and real term of electrostatic interaction
  overlap = calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
  if(!overlap) {
    //calculate reciprocate term of electrostatic interaction
    recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
  }
}

inline void Rotate::Accept(const uint rejectState, const uint step)
{
  bool res = false;

  if(rejectState == mv::fail_state::NO_FAIL) {
    real pr = prng();
    res = pr < exp(-BETA * (inter_LJ.energy + inter_Real.energy +
                            recip.energy));
  }
  bool result = res && !overlap;

  if (result) {
    //Set new energy.
    // setting energy and virial of LJ interaction
    sysPotRef.boxEnergy[b].inter += inter_LJ.energy;
    // setting energy and virial of coulomb interaction
    sysPotRef.boxEnergy[b].real += inter_Real.energy;
    // setting energy and virial of recip term
    sysPotRef.boxEnergy[b].recip += recip.energy;

    //Copy coords
    newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
    calcEwald->UpdateRecip(b);

    sysPotRef.Total();
  }

  if(molRemoved) {
    // It means that Recip energy is calculated and move not accepted
    if(!result && !overlap) {
      calcEwald->RestoreMol(m);
    }

    cellList.AddMol(m, b, coordCurrRef);
    molRemoved = false;
  }

  moveSetRef.Update(mv::ROTATE, result, step, b, mk);
}

#endif /*ROTATION_H*/
