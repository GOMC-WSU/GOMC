#pragma once

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"

#define MPDISPLACE 0
#define MPROTATE 1

class MultiParticle : public MoveBase
{
public:
  MultiParticle(System &sys, StaticVals const& statV);

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  double GetCoeff();
private:
  uint bPick;
  double w_new, w_old;
  double t_max, r_max;
  double lambda;
  SystemPotential sysPotNew;
  XYZArray& atomTorqueRef;
  XYZArray atomTorqueNew;
  XYZArray& molTorqueRef;
  XYZArray molTorqueNew;
  XYZArray t_k;
  XYZArray r_k;
  Coordinates newMolsPos;
  COM newCOMs;
  vector<uint> moveType;
  const MoleculeLookup& molLookup;

  void CalculateTrialDistRot(uint molIndex);
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV),
  newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef,statV.mol),
  atomTorqueRef(sys.atomTorqueRef),
  molTorqueRef(sys.molTorqueRef),
  molLookup(sys.molLookup)
{
  atomTorqueNew.Init(sys.atomTorqueRef.Count());
  molTorqueNew.Init(sys.molTorqueRef.Count());
  t_k.Init(sys.com.Count());
  r_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());
  moveType.resize(sys.com.Count());
  lambda = 0.5;
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif
  // subPick = mv::GetMoveSubIndex(mv::MULTIPARTICLE, bPick);
  t_max = 0.05;
  r_max = 0.09 * 2 * M_PI;

  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(bPick);
  while(thisMol != end) {
    uint length = molRef.GetKind(*thisMol).NumAtoms();
    if(length==1)
      moveType[*thisMol] = MPDISPLACE;
    else
      moveType[*thisMol] = prng.randInt(1);
    CalculateTrialDistRot(*thisMol);
    thisMol++;
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  return state;
}

inline uint MultiParticle::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.

  // move particles according to force and torque and store them in the new pos
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(bPick);
  while(thisMol != end) {
    double molIndex = (*thisMol);
    if(moveType[molIndex]) { // rotate
      RotateForceBiased(molIndex);
    } else { // displacement
      TranslateForceBiased(molIndex);
    }
    thisMol++;
  } 
  return 0;
}

inline void MultiParticle::CalcEn() 
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  cellList.GridAll(boxDimRef, newMolsPos, molLookup);

  sysPotNew = sysPotRef;
  sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newCOMs, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  calcEnRef.CalculateTorque(coordCurrRef, comCurrRef, atomForceRef, atomTorqueRef,
                            molTorqueRef, bPick);
  calcEnRef.CalculateTorque(newMolsPos, newCOMs, atomForceNew, atomTorqueNew,
                            molTorqueNew, bPick);
  return;
}

inline double MultiParticle::GetCoeff()
{
  // calculate (w_new->old/w_old->new) and return it.
  uint length, start;
  XYZ lbf; // lambda * BETA * force
  XYZ lbt; // lambda * BETA * torque
  w_old = 1.0;
  w_new = 1.0;
  uint molNumber;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(bPick);

  while(thisMol != end) {
    molNumber = *thisMol;
    if(moveType[molNumber]) { // == 1 -> rotate
      lbt = molTorqueRef.Get(molNumber) * lambda * BETA;
      w_new *= lbt.x * exp(lbf.x * r_k.Get(molNumber).x)/
        (2.0*sinh(lbt.x * r_max));
      w_new *= lbt.y * exp(lbf.y * r_k.Get(molNumber).y)/
        (2.0*sinh(lbt.y * r_max));
      w_new *= lbt.z * exp(lbf.z * r_k.Get(molNumber).z)/
        (2.0*sinh(lbt.z * r_max));

      lbt = molTorqueNew.Get(molNumber) * lambda * BETA;
      w_old *= lbt.x * exp(lbf.x * -1 * r_k.Get(molNumber).x)/
        (2.0*sinh(lbt.x * r_max));
      w_old *= lbt.y * exp(lbf.y * -1 * r_k.Get(molNumber).y)/
        (2.0*sinh(lbt.y * r_max));
      w_old *= lbt.z * exp(lbf.z * -1 * r_k.Get(molNumber).z)/
        (2.0*sinh(lbt.z * r_max));
    }
    else { // displace
      lbf = molForceRef.Get(molNumber) * lambda * BETA;
      w_new *= lbf.x * exp(lbf.x * t_k.Get(molNumber).x)/
        (2.0*sinh(lbf.x * t_max));
      w_new *= lbf.y * exp(lbf.y * t_k.Get(molNumber).y)/
        (2.0*sinh(lbf.y * t_max));
      w_new *= lbf.z * exp(lbf.z * t_k.Get(molNumber).z)/
        (2.0*sinh(lbf.z * t_max));

      lbf = molForceNew.Get(molNumber) * lambda * BETA;
      w_old *= lbf.x * exp(lbf.x * -1 * t_k.Get(molNumber).x)/
        (2.0*sinh(lbf.x * t_max));
      w_old *= lbf.y * exp(lbf.y * -1 * t_k.Get(molNumber).y)/
        (2.0*sinh(lbf.y * t_max));
      w_old *= lbf.z * exp(lbf.z * -1 * t_k.Get(molNumber).z)/
        (2.0*sinh(lbf.z * t_max));

    }
    thisMol++;
  }
  return w_new/w_old;
}

inline void MultiParticle::Accept(const uint rejectState, const uint step)
{
  // Here we compare the values of reference and trial and decide whether to 
  // accept or reject the move
  double MPCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = MPCoeff * uBoltz;
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if(result) {
    sysPotRef = sysPotNew;
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    swap(molForceRef, molForceNew);
    swap(atomForceRef, atomForceNew);
    swap(molTorqueRef, molTorqueNew);
    swap(atomTorqueRef, atomTorqueNew);
  }
  else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookup);
  }
  subPick = mv::GetMoveSubIndex(mv::MULTIPARTICLE, bPick);
  moveSetRef.Update(result, subPick, step);
}

inline void MultiParticle::CalculateTrialDistRot(uint molIndex)
{
  XYZ lbf, lbfmax; // lambda * BETA * force
  XYZ lbt, lbtmax; // lambda * BETA * torque
  double rand;
  XYZ num;
  if(moveType[molIndex]) { // rotate
    lbt = molTorqueRef.Get(molIndex) * lambda * BETA;
    lbtmax = lbt * r_max;
    rand = prng();
    num.x = log(exp(-1 * lbtmax.x ) + 2 * rand * sinh(lbtmax.x ));
    rand = prng();
    num.y = log(exp(-1 * lbtmax.y ) + 2 * rand * sinh(lbtmax.y ));
    rand = prng();
    num.z = log(exp(-1 * lbtmax.z ) + 2 * rand * sinh(lbtmax.z ));
    r_k.Set(molIndex, num/lbt);
  }
  else { // displace
    lbf = molForceRef.Get(molIndex) * lambda * BETA;
    lbfmax = lbf * t_max;
    rand = prng();
    num.x = log(exp(-1 * lbfmax.x ) + 2 * rand * sinh(lbfmax.x ));
    rand = prng();
    num.y = log(exp(-1 * lbfmax.y ) + 2 * rand * sinh(lbfmax.y ));
    rand = prng();
    num.z = log(exp(-1 * lbfmax.z ) + 2 * rand * sinh(lbfmax.z ));
    t_k.Set(molIndex, num/lbf);
  }
}

void MultiParticle::RotateForceBiased(uint molIndex)
{
  
}

void MultiParticle::TranslateForceBiased(uint molIndex)
{
  XYZ shift = t_k.Get(molIndex);
  uint stop, start, len;

  molRef.GetRange(start, stop, len, molIndex);
  //Add translation
  newMolsPos.AddRange(start, stop, shift);
  newCOMs.Set(molIndex, newCOMs.Get(molIndex) + shift);

  for(uint i=start; i<stop; i++)
    newMolsPos.Set(i, boxDimRef.WrapPBC(newMolsPos.Get(i), bPick));
  newCOMs.Set(molIndex, boxDimRef.WrapPBC(newMolsPos.Get(molIndex), bPick));
}