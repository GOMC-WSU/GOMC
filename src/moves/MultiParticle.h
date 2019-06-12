/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MULTIPARTICLE_H
#define MULTIPARTICLE_H

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"
#include <cmath>

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

class MultiParticle : public MoveBase
{
public:
  MultiParticle(System &sys, StaticVals const& statV);

  virtual uint Prep(const real subDraw, const real movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();

private:
  uint bPick;
  uint typePick;
  real lambda;
  bool initMol[BOX_TOTAL];
  SystemPotential sysPotNew;
  XYZArray molTorqueRef;
  XYZArray molTorqueNew;
  XYZArray atomForceRecNew;
  XYZArray molForceRecNew;
  XYZArray t_k;
  XYZArray r_k;
  Coordinates newMolsPos;
  COM newCOMs;
  vector<uint> moveType, moleculeIndex;
  const MoleculeLookup& molLookup;

  long real GetCoeff();
  void CalculateTrialDistRot();
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(XYZ const &lb, real const max);
  real CalculateWRatio(XYZ const &lb, XYZ const &k, real max,
                            real sign);
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

  // set default value for r_max, t_max, and lambda
  // the value of lambda is based on the paper
  lambda = 0.5;
  for(uint b = 0; b < BOX_TOTAL; b++) {
    initMol[b] = false;
  }
}

inline void MultiParticle::PrintAcceptKind() {
  printf("%-37s", "% Accepted MultiParticle ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::MULTIPARTICLE));
  }
  std::cout << std::endl;
}


inline void MultiParticle::SetMolInBox(uint box)
{
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  moleculeIndex.clear();
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
  while(thisMol != end) {
    moleculeIndex.push_back(*thisMol);
    thisMol++;
  }
#else
  if(!initMol[box]){  
    moleculeIndex.clear();
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    while(thisMol != end) {
      moleculeIndex.push_back(*thisMol);
      thisMol++;
    }
  }
#endif
  initMol[box] = true;
}

inline uint MultiParticle::Prep(const real subDraw, const real movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif

  typePick = prng.randIntExc(mp::MPTOTALTYPES);
  SetMolInBox(bPick);

  for(uint m = 0; m < moleculeIndex.size(); m++) {
    uint length = molRef.GetKind(moleculeIndex[m]).NumAtoms();
    if(length == 1) {
      moveType[moleculeIndex[m]] = mp::MPDISPLACE;
    } else {
      switch(typePick) {
        case mp::MPALLDISPLACE:
          moveType[moleculeIndex[m]] = mp::MPDISPLACE;
          break;
        case mp::MPALLROTATE:
          moveType[moleculeIndex[m]] = mp::MPROTATE;
          break;
        case mp::MPALLRANDOM:
          moveType[moleculeIndex[m]] = (prng.randInt(1) ?
					mp::MPROTATE : mp::MPDISPLACE);
          break;
        default:
          std::cerr << "Error: Something went wrong preping MultiParticle!\n"
                    << "Type Pick value is not valid!\n";
          exit(EXIT_FAILURE);
      }
    }
  }

  if(moveSetRef.GetSingleMoveAccepted()){
    //Calculate force for long range electrostatic using old position
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef,molForceRecRef,
                                  bPick);
    
    //calculate short range energy and force for old positions
    calcEnRef.BoxInter(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                      boxDimRef, bPick);

    if(typePick != mp::MPALLDISPLACE) {
      //Calculate Torque for old positions
      calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                                atomForceRef, atomForceRecRef, molTorqueRef,
                                moveType, bPick);
    }
  }
  CalculateTrialDistRot();
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  return state;
}

inline uint MultiParticle::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;
  uint m;

  // move particles according to force and torque and store them in the new pos
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(m)
#endif
  for(m = 0; m < moleculeIndex.size(); m++) {
    if(moveType[moleculeIndex[m]]) {
      // rotate
      RotateForceBiased(moleculeIndex[m]);
    } else {
      // displacement
      TranslateForceBiased(moleculeIndex[m]);
    }
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
  sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  //calculate long range of new electrostatic energy
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick);
  //Calculate long range of new electrostatic force
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecNew, molForceRecNew,
    bPick);
  //Calculate Torque for new positions
  calcEnRef.CalculateTorque(moleculeIndex, newMolsPos, newCOMs, atomForceNew,
                            atomForceRecNew, molTorqueNew, moveType, bPick);
  sysPotNew.Total();
}

inline real MultiParticle::CalculateWRatio(XYZ const &lb, XYZ const &k,
                                                real max, real sign)
{
  real w_ratio = 1.0;
  XYZ lbmax = lb * max;

  if(abs(lbmax.x) > MIN_FORCE && abs(lbmax.x) < MAX_FORCE) {
    w_ratio *= lb.x * exp(lb.x * sign * k.x) / (2.0*sinh(lb.x * max));
  }  else {
    w_ratio *= 1.0 / (2.0 * max);
  }
  
  if(abs(lbmax.y) > MIN_FORCE && abs(lbmax.y) < MAX_FORCE){
    w_ratio *= lb.y * exp(lb.y * sign * k.y) / (2.0*sinh(lb.y * max));
  } else {
    w_ratio *= 1.0 / (2.0 * max);
  }

  if(abs(lbmax.z) > MIN_FORCE && abs(lbmax.z) < MAX_FORCE){
    w_ratio *= lb.z * exp(lb.z * sign * k.z) / (2.0*sinh(lb.z * max));
  } else {
    w_ratio *= 1.0 / (2.0 * max);
  }

  // if(sign == -1) {
  //   cout << "w_ratio: " << w_ratio << endl;
  // }
  return w_ratio;
}

inline long real MultiParticle::GetCoeff()
{
  // calculate (w_new->old/w_old->new) and return it.
  XYZ lbf_old, lbf_new; // lambda * BETA * force
  XYZ lbt_old, lbt_new; // lambda * BETA * torque
  long real w_ratio = 1.0;
  real lBeta = lambda * BETA;
  uint m, molNumber;
  real r_max = moveSetRef.GetRMAX(bPick);
  real t_max = moveSetRef.GetTMAX(bPick);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(m, molNumber, lbt_old, lbt_new, lbf_old, lbf_new) reduction(*:w_ratio)
#endif
  for(m = 0; m < moleculeIndex.size(); m++) {
    molNumber = moleculeIndex[m];
    if(moveType[molNumber]) {
      // rotate
      lbt_old = molTorqueRef.Get(molNumber) * lBeta;
      lbt_new = molTorqueNew.Get(molNumber) * lBeta;

      w_ratio *= CalculateWRatio(lbt_new, r_k.Get(molNumber), r_max, -1);
      w_ratio /= CalculateWRatio(lbt_old, r_k.Get(molNumber), r_max, 1);
    } else {
      // displace
      lbf_old = (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) *
	      lBeta;
      lbf_new = (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) *
        lBeta;
      w_ratio *= CalculateWRatio(lbf_new, t_k.Get(molNumber), t_max, -1);
      w_ratio /= CalculateWRatio(lbf_old, t_k.Get(molNumber), t_max, 1);
    }
  }

  // In case where force or torque is a large negative number (ex. -800)
  // the exp value becomes inf. In these situations we have to return 0 to
  // reject the move
  // if(!std::isfinite(w_ratio)) {
  //   // This error can be removed later on once we know this part actually works.
  //   std::cout << "w_ratio is not a finite number. Auto-rejecting move.\n";
  //   return 0.0;
  // }
  return w_ratio;
}

inline void MultiParticle::Accept(const uint rejectState, const uint step)
{
  // Here we compare the values of reference and trial and decide whether to 
  // accept or reject the move
  long real MPCoeff = GetCoeff();
  real uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  long real accept = MPCoeff * uBoltz;
  // cout << "MPCoeff: " << MPCoeff << ", sysPotNew: " << sysPotNew.Total()
  //      << ", sysPotRef: " << sysPotRef.Total() << ", accept: " << accept <<endl;
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

  moveSetRef.UpdateMoveSettingMultiParticle(bPick, result, typePick);
  moveSetRef.AdjustMultiParticle(bPick, typePick);

  moveSetRef.Update(mv::MULTIPARTICLE, result, step, bPick);
}

inline XYZ MultiParticle::CalcRandomTransform(XYZ const &lb, real const max)
{
  XYZ lbmax = lb * max;
  XYZ num;
  if(abs(lbmax.x) > MIN_FORCE && abs(lbmax.x) < MAX_FORCE) {
    num.x = log(exp(-1.0 * lbmax.x ) + 2 * prng() * sinh(lbmax.x )) / lb.x;
  } else {
    num.x = max * prng.Sym(1);
  }
  
  if(abs(lbmax.y) > MIN_FORCE && abs(lbmax.y) < MAX_FORCE){
    num.y = log(exp(-1.0 * lbmax.y ) + 2 * prng() * sinh(lbmax.y )) / lb.y;
  } else {
    num.y = max * prng.Sym(1);
  }

  if(abs(lbmax.z) > MIN_FORCE && abs(lbmax.z) < MAX_FORCE){
    num.z = log(exp(-1.0 * lbmax.z ) + 2 * prng() * sinh(lbmax.z )) / lb.z;
  } else {
    num.z = max * prng.Sym(1);
  }

  if(num.Length() >= boxDimRef.axis.Min(bPick)) {
    std::cout << "Trial Displacement exceed half of the box length in Multiparticle move.\n";
    std::cout << "Trial transform: " << num;
    exit(EXIT_FAILURE);
  } else if (!isfinite(num.Length())) {
    std::cout << "Trial Displacement is not a finite number in Multiparticle move.\n";
    std::cout << "Trial transform: " << num;
    exit(EXIT_FAILURE);
  }

  // We can possible bound them

  return num;
}

inline void MultiParticle::CalculateTrialDistRot()
{
  uint m , molIndex;
  real r_max = moveSetRef.GetRMAX(bPick);
  real t_max = moveSetRef.GetTMAX(bPick);
  XYZ lbf; // lambda * BETA * force * maxTranslate
  XYZ lbt; // lambda * BETA * torque * maxRotation
  for(m = 0; m < moleculeIndex.size(); m++) {
    molIndex = moleculeIndex[m];

    if(moveType[molIndex]) { // rotate
      lbt = molTorqueRef.Get(molIndex) * lambda * BETA;
      r_k.Set(molIndex, CalcRandomTransform(lbt, r_max));
    } else { // displace
      lbf = (molForceRef.Get(molIndex) + molForceRecRef.Get(molIndex)) *
        lambda * BETA;
      t_k.Set(molIndex, CalcRandomTransform(lbf, t_max));
    }
  }
}

inline void MultiParticle::RotateForceBiased(uint molIndex)
{
  XYZ rot = r_k.Get(molIndex);
  real rotLen = rot.Length();
  RotationMatrix matrix;
  
  matrix = RotationMatrix::FromAxisAngle(rotLen, rot * (1.0/rotLen));
 
  XYZ center = newCOMs.Get(molIndex);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);
  
  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);
  
  // Do Rotation
  for(uint p=0; p<len; p++) {
    temp.Add(p, -center);
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, center);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticle::TranslateForceBiased(uint molIndex)
{
  XYZ shift = t_k.Get(molIndex);

  XYZ newcom = newCOMs.Get(molIndex);
  uint stop, start, len;
  molRef.GetRange(start, stop, len, molIndex);
  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  //Shift the coordinate and COM
  temp.AddAll(shift);
  newcom += shift;
  //rewrapping
  boxDimRef.WrapPBC(temp, bPick);
  newcom = boxDimRef.WrapPBC(newcom, bPick);
  //set the new coordinate
  temp.CopyRange(newMolsPos, 0, start, len);
  newCOMs.Set(molIndex, newcom);
}

#endif
