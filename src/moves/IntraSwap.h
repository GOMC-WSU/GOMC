/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef INTRASWAP_H
#define INTRASWAP_H


#include "MoveBase.h"
#include "TrialMol.h"

//#define DEBUG_MOVES

class IntraSwap : public MoveBase
{
public:

  IntraSwap(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
    MoveBase(sys, statV) {}

  virtual uint Prep(const real subDraw, const real movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);
  virtual void PrintAcceptKind();
private:
  uint GetBoxAndMol(const real subDraw, const real movPerc);
  MolPick molPick;
  uint sourceBox, destBox;
  uint pStart, pLen;
  uint molIndex, kindIndex;

  real W_tc, W_recip;
  real correct_old, correct_new;
  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipDiff;
  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;
};

void IntraSwap::PrintAcceptKind()
{
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted Intra-Swap ", molRef.kinds[k].name.c_str());
    for(uint b = 0; b < BOX_TOTAL; b++) {
      if(moveSetRef.GetTrial(b, mv::INTRA_SWAP, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::INTRA_SWAP, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint IntraSwap::GetBoxAndMol(const real subDraw, const real movPerc)
{

#if ENSEMBLE == GCMC
  sourceBox = mv::BOX0;
  uint state = prng.PickMol(molIndex, kindIndex, sourceBox, subDraw, movPerc);
#else
  uint state = prng.PickMolAndBox(molIndex, kindIndex, sourceBox, subDraw,
                                  movPerc);
#endif

  //molecule will be removed and insert in same box
  destBox = sourceBox;

  if (state == mv::fail_state::NO_FAIL) {
    pStart = pLen = 0;
    molRef.GetRangeStartLength(pStart, pLen, molIndex);
  }
  return state;
}

inline uint IntraSwap::Prep(const real subDraw, const real movPerc)
{
  overlap = false;
  uint state = GetBoxAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStart);
    W_tc = 1.0;
  }
  return state;
}


inline uint IntraSwap::Transform()
{
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);
  overlap = newMol.HasOverlap();
  return mv::fail_state::NO_FAIL;
}

inline void IntraSwap::CalcEn()
{
  // since number of molecules would not change in the box,
  //there is no change in Tc
  W_tc = 1.0;
  W_recip = 1.0;
  correct_old = 0.0;
  correct_new = 0.0;

  if (newMol.GetWeight() != 0.0 && !overlap) {
    correct_new = calcEwald->SwapCorrection(newMol);
    correct_old = calcEwald->SwapCorrection(oldMol);
    recipDiff.energy = calcEwald->MolReciprocal(newMol.GetCoords(), molIndex,
                       sourceBox);
    //self energy is same
    W_recip = exp(-1.0 * ffRef.beta * (recipDiff.energy + correct_new -
                                       correct_old));
  }
}


inline void IntraSwap::Accept(const uint rejectState, const uint step)
{
  bool result;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    real molTransCoeff = 1.0;
    real Wo = oldMol.GetWeight();
    real Wn = newMol.GetWeight();
    real Wrat = Wn / Wo * W_tc * W_recip;

    //safety to make sure move will be rejected in overlap case
    if(!overlap) {
      result = prng() < molTransCoeff * Wrat;
    } else
      result = false;

    if(result) {
      //Add rest of energy.
      sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
      sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
      //Add Reciprocal energy difference
      sysPotRef.boxEnergy[destBox].recip += recipDiff.energy;
      //Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
      sysPotRef.boxEnergy[destBox].correction += correct_new;

      //Set coordinates, new COM; shift index to new box's list
      newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
      comCurrRef.SetNew(molIndex, destBox);
      cellList.AddMol(molIndex, destBox, coordCurrRef);

      //Zero out box energies to prevent small number
      //errors in double.
      if (molLookRef.NumInBox(sourceBox) == 1) {
        sysPotRef.boxEnergy[sourceBox].inter = 0;
        sysPotRef.boxVirial[sourceBox].inter = 0;
        sysPotRef.boxEnergy[sourceBox].real_en = 0;
        sysPotRef.boxVirial[sourceBox].real_en = 0;
      }

      calcEwald->UpdateRecip(destBox);
      
      //Retotal
      sysPotRef.Total();
    } else {
      cellList.AddMol(molIndex, sourceBox, coordCurrRef);

      //when weight is 0, MolDestSwap() will not be executed, thus cos/sin
      //molRef will not be changed. Also since no memcpy, doing restore
      //results in memory overwrite
      if (newMol.GetWeight() != 0.0 && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }
  } else //else we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::INTRA_SWAP, result, step, sourceBox, kindIndex);
}

#endif
