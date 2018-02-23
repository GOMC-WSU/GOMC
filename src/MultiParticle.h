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
  XYZArray molTorqueRef;
  XYZArray molTorqueNew;
  XYZArray atomForceRecNew;
  XYZArray molForceRecNew;
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
  molLookup(sys.molLookup)
{
  molTorqueNew.Init(sys.com.Count());
  molTorqueRef.Init(sys.com.Count());
  atomForceRecNew.Init(sys.coordinates.Count());
  molForceRecNew.Init(sys.com.Count());
  
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
  uint state = mv::fail_state::NO_FAIL;

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
  return state;
}

inline void MultiParticle::CalcEn() 
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  cellList.GridAll(boxDimRef, newMolsPos, molLookup);

  //back up cached fourier term
  calcEwald->exgMolCache();
  //setup reciprocate vectors for new positions
  calcEwald->BoxReciprocalSetup(bPick, newMolsPos);

  sysPotNew = sysPotRef;
  //calculate short range energy and force
  sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newCOMs, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  //calculate long range of new electrostatic energy
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick);
  //Calculate long range of new electrostatic force
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecNew, molForceRecNew,
				bPick);
  //Calculate Torque for old positions
  calcEnRef.CalculateTorque(coordCurrRef, comCurrRef, atomForceRef,
			    atomForceRecRef, molTorqueRef, moveType, bPick);
  //Calculate Torque for new positions
  calcEnRef.CalculateTorque(newMolsPos, newCOMs, atomForceNew, atomForceRecNew,
			    molTorqueNew, moveType, bPick);
  sysPotNew.Total();
}

inline double MultiParticle::GetCoeff()
{
  // calculate (w_new->old/w_old->new) and return it.
  uint length, start;
  XYZ lbf_old, lbf_new; // lambda * BETA * force
  XYZ lbt_old, lbt_new; // lambda * BETA * torque
  double w_ratio = 1.0;
  uint molNumber;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(bPick);

  while(thisMol != end) {
    molNumber = *thisMol;
    if(moveType[molNumber]) { // == 1 -> rotate
      lbt_old = molTorqueRef.Get(molNumber) * lambda * BETA;
      lbt_new = molTorqueNew.Get(molNumber) * lambda * BETA;
      w_ratio *= lbt_new.x * exp(lbt_new.x * -1 * r_k.Get(molNumber).x)/
        (2.0*sinh(lbt_new.x * r_max));
      w_ratio *= lbt_new.y * exp(lbt_new.y * -1 * r_k.Get(molNumber).y)/
        (2.0*sinh(lbt_new.y * r_max));
      w_ratio *= lbt_new.z * exp(lbt_new.z * -1 * r_k.Get(molNumber).z)/
        (2.0*sinh(lbt_new.z * r_max));

      w_ratio /= lbt_old.x * exp(lbt_old.x * r_k.Get(molNumber).x)/
        (2.0*sinh(lbt_old.x * r_max));
      w_ratio /= lbt_old.y * exp(lbt_old.y * r_k.Get(molNumber).y)/
        (2.0*sinh(lbt_old.y * r_max));
      w_ratio /= lbt_old.z * exp(lbt_old.z * r_k.Get(molNumber).z)/
        (2.0*sinh(lbt_old.z * r_max));
    }
    else { // displace
      lbf_old = (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) *
	lambda * BETA;
      lbf_new = (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) *
	lambda * BETA;
      w_ratio *= lbf_new.x * exp(lbf_new.x * -1 * t_k.Get(molNumber).x)/
        (2.0*sinh(lbf_new.x * t_max));
      w_ratio *= lbf_new.y * exp(lbf_new.y * -1 * t_k.Get(molNumber).y)/
        (2.0*sinh(lbf_new.y * t_max));
      w_ratio *= lbf_new.z * exp(lbf_new.z * -1 * t_k.Get(molNumber).z)/
        (2.0*sinh(lbf_new.z * t_max));

      w_ratio /= lbf_old.x * exp(lbf_old.x * t_k.Get(molNumber).x)/
        (2.0*sinh(lbf_old.x * t_max));
      w_ratio /= lbf_old.y * exp(lbf_old.y * t_k.Get(molNumber).y)/
        (2.0*sinh(lbf_old.y * t_max));
      w_ratio /= lbf_old.z * exp(lbf_old.z * t_k.Get(molNumber).z)/
        (2.0*sinh(lbf_old.z * t_max));
    }
    thisMol++;
  }
  return w_ratio;
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
    swap(molForceRecRef, molForceRecNew);
    swap(atomForceRecRef, atomForceRecNew);
    swap(molTorqueRef, molTorqueNew);
    //update reciprocate value
    calcEwald->UpdateRecip(bPick);
  }
  else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookup);
    calcEwald->exgMolCache();
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
    XYZ temp = num/lbt;
    if(isnan(temp.x))
      temp.x = 0.0;
    if(isnan(temp.y))
      temp.y = 0.0;
    if(isnan(temp.z))
      temp.z = 0.0;
    r_k.Set(molIndex, temp);
  }
  else { // displace
    lbf = (molForceRef.Get(molIndex) * molForceRecRef.Get(molIndex)) *
      lambda * BETA;
    lbfmax = lbf * t_max;
    rand = prng();
    num.x = log(exp(-1 * lbfmax.x ) + 2 * rand * sinh(lbfmax.x ));
    rand = prng();
    num.y = log(exp(-1 * lbfmax.y ) + 2 * rand * sinh(lbfmax.y ));
    rand = prng();
    num.z = log(exp(-1 * lbfmax.z ) + 2 * rand * sinh(lbfmax.z ));
    XYZ temp = num/lbf;
    if(isnan(temp.x))
      temp.x = 0.0;
    if(isnan(temp.y))
      temp.y = 0.0;
    if(isnan(temp.z))
      temp.z = 0.0;
    t_k.Set(molIndex, temp);
  }
}

void MultiParticle::RotateForceBiased(uint molIndex)
{
  XYZ rot = r_k.Get(molIndex);
  RotationMatrix matrix = RotationMatrix::FromAxisAngle(rot.Length(), rot);
  XYZ center = newCOMs.Get(molIndex);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  for(uint i=start; i<stop; i++) {
    XYZ ith = newMolsPos[i];
    XYZ unwrapped = boxDimRef.UnwrapPBC(ith, bPick, center);
    newMolsPos.Set(i, unwrapped);
    newMolsPos.Add(i, -center);
    newMolsPos.Set(i, matrix.Apply(newMolsPos[i]));
    newMolsPos.Add(i, center);
    newMolsPos.Set(i, boxDimRef.WrapPBC(newMolsPos[i], bPick));
  }
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
