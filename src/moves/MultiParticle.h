/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
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
#include <fstream>
#include "Random123Wrapper.h"
#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

class MultiParticle : public MoveBase
{
public:
  MultiParticle(System &sys, StaticVals const& statV);

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();
  void PrepCFCMC(const uint box);

private:
  uint bPick;
  double lambda;
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
  int moveType;
  std::vector<uint> moleculeIndex;
  const MoleculeLookup& molLookup;
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars;
  std::vector<int> particleMol;
#endif
  Random123Wrapper &r123wrapper;
  const Molecules& mols;

  double GetCoeff();
  void CalculateTrialDistRot();
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(XYZ const &lb, double const max, uint molIndex);
  double CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old, XYZ const &k,
                         double max);
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV),

  newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol),
  molLookup(sys.molLookup), r123wrapper(sys.r123wrapper), mols(statV.mol)
{
  molTorqueNew.Init(sys.com.Count());
  molTorqueRef.Init(sys.com.Count());
  atomForceRecNew.Init(sys.coordinates.Count());
  molForceRecNew.Init(sys.com.Count());

  t_k.Init(sys.com.Count());
  r_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());

  // set default value for r_max, t_max, and lambda
  // the value of lambda is based on the paper
  lambda = 0.5;
  for(uint b = 0; b < BOX_TOTAL; b++) {
    initMol[b] = false;
  }

#ifdef GOMC_CUDA
  cudaVars = sys.statV.forcefield.particles->getCUDAVars();

  uint maxAtomInMol = 0;
  for(uint m = 0; m < mols.count; ++m) {
    const MoleculeKind& molKind = mols.GetKind(m);
    if(molKind.NumAtoms() > maxAtomInMol)
      maxAtomInMol = molKind.NumAtoms();
    for(uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleMol.push_back(m);
    }
  }
#endif
}

inline void MultiParticle::PrintAcceptKind()
{
  printf("%-37s", "% Accepted MultiParticle ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::MULTIPARTICLE));
  }
  std::cout << std::endl;
}


inline void MultiParticle::SetMolInBox(uint box)
{
  // NEED to check if atom is not fixed!
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  moleculeIndex.clear();
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
  while(thisMol != end) {
    //Make sure this molecule is not fixed in its position
    if(!molLookup.IsFix(*thisMol)) {
      moleculeIndex.push_back(*thisMol);
    }
    thisMol++;
  }
#else
  if(!initMol[box]) {
    moleculeIndex.clear();
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    while(thisMol != end) {
      //Make sure this molecule is not fixed in its position
      if(!molLookup.IsFix(*thisMol)) {
        moleculeIndex.push_back(*thisMol);
      }
      thisMol++;
    }
  }
#endif
  initMol[box] = true;
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif

  moveType = prng.randIntExc(mp::MPTOTALTYPES);
  SetMolInBox(bPick);
  if (moleculeIndex.size() == 0) {
    std::cout << "Warning: Multi particle move can't move any molecules, Skipping...\n";
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    return state;
  }
  uint length = molRef.GetKind(moleculeIndex[0]).NumAtoms();
  if(length == 1) {
    moveType = mp::MPDISPLACE;
  }

  if(moveSetRef.GetSingleMoveAccepted()) {
    //Calculate force for long range electrostatic using old position
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                  bPick);
    //calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);
    if(moveType == mp::MPROTATE) {
      //Calculate Torque for old positions
      calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                                atomForceRef, atomForceRecRef, molTorqueRef,
                                bPick);
    }
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  return state;
}

inline void MultiParticle::PrepCFCMC(const uint box)
{
  bPick = box;
  moveType = prng.randIntExc(mp::MPTOTALTYPES);
  SetMolInBox(bPick);
  uint length = molRef.GetKind(moleculeIndex[0]).NumAtoms();
  if(length == 1) {
    moveType = mp::MPDISPLACE;
  }

  if(moveSetRef.GetSingleMoveAccepted()) {
    //Calculate force for long range electrostatic using old position
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                  bPick);

    //calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);

    if(moveType == mp::MPROTATE) {
      //Calculate Torque for old positions
      calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                                atomForceRef, atomForceRecRef, molTorqueRef,
                                bPick);
    }
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
}

inline uint MultiParticle::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;
#ifdef GOMC_CUDA
  // calculate which particles are inside moleculeIndex
  std::vector<int> isMoleculeInvolved(newMolsPos.Count(), 0);
  for(int particleID = 0; particleID < isMoleculeInvolved.size(); particleID++) {
    int midx = particleMol[particleID];
    std::vector<uint>::iterator it;
    it = find(moleculeIndex.begin(), moleculeIndex.end(), midx);
    if(it != moleculeIndex.end()) {
      isMoleculeInvolved[particleID] = 1;
    }
  }

  // This kernel will calculate translation/rotation amount + shifting/rotating
  if(moveType == mp::MPROTATE) {
    double r_max = moveSetRef.GetRMAX(bPick);
    CallRotateParticlesGPU(cudaVars, isMoleculeInvolved, r_max,
                           molTorqueRef.x, molTorqueRef.y, molTorqueRef.z,
                           r123wrapper.GetStep(), r123wrapper.GetSeedValue(),
                           particleMol, atomForceRecNew.Count(),
                           molForceRecNew.Count(), boxDimRef.GetAxis(bPick).x,
                           boxDimRef.GetAxis(bPick).y, boxDimRef.GetAxis(bPick).z,
                           newMolsPos, newCOMs, lambda * BETA, r_k);
  } else {
    double t_max = moveSetRef.GetTMAX(bPick);
    CallTranslateParticlesGPU(cudaVars, isMoleculeInvolved, t_max,
                              molForceRef.x, molForceRef.y, molForceRef.z,
                              r123wrapper.GetStep(), r123wrapper.GetSeedValue(),
                              particleMol, atomForceRecNew.Count(),
                              molForceRecNew.Count(), boxDimRef.GetAxis(bPick).x,
                              boxDimRef.GetAxis(bPick).y, boxDimRef.GetAxis(bPick).z,
                              newMolsPos, newCOMs, lambda * BETA, t_k, molForceRecRef);
  }
#else
  CalculateTrialDistRot();
  // move particles according to force and torque and store them in the new pos
  if(moveType == mp::MPROTATE) {
#ifdef _OPENMP
    #pragma omp parallel for default(none)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      RotateForceBiased(moleculeIndex[m]);
    }
  } else {
#ifdef _OPENMP
    #pragma omp parallel for default(none)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      TranslateForceBiased(moleculeIndex[m]);
    }
  }
#endif
  return state;
}

inline void MultiParticle::CalcEn()
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  cellList.GridAll(boxDimRef, newMolsPos, molLookup);

  //back up cached Fourier term
  calcEwald->backupMolCache();
  //setup reciprocal vectors for new positions
  calcEwald->BoxReciprocalSetup(bPick, newMolsPos);

  sysPotNew = sysPotRef;
  //calculate short range energy and force
  sysPotNew = calcEnRef.BoxForce(sysPotNew, newMolsPos, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  //calculate long range of new electrostatic energy
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick);
  //Calculate long range of new electrostatic force
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecNew, molForceRecNew,
                                bPick);

  if(moveType == mp::MPROTATE) {
    //Calculate Torque for new positions
    calcEnRef.CalculateTorque(moleculeIndex, newMolsPos, newCOMs, atomForceNew,
                              atomForceRecNew, molTorqueNew, bPick);
  }
  sysPotNew.Total();
}

inline double MultiParticle::CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old,
    XYZ const &k, double max)
{
  double w_ratio = 1.0;
  XYZ lbmax = lb_old * max;
  //If we used force to bias the displacement or rotation, we include it
  if(std::abs(lbmax.x) > MIN_FORCE && std::abs(lbmax.x) < MAX_FORCE) {
    w_ratio *= lb_new.x * exp(-lb_new.x * k.x) / (2.0 * sinh(lb_new.x * max));
    w_ratio /= lb_old.x * exp(lb_old.x * k.x) / (2.0 * sinh(lb_old.x * max));
  }

  if(std::abs(lbmax.y) > MIN_FORCE && std::abs(lbmax.y) < MAX_FORCE) {
    w_ratio *= lb_new.y * exp(-lb_new.y * k.y) / (2.0 * sinh(lb_new.y * max));
    w_ratio /= lb_old.y * exp(lb_old.y * k.y) / (2.0 * sinh(lb_old.y * max));
  }

  if(std::abs(lbmax.z) > MIN_FORCE && std::abs(lbmax.z) < MAX_FORCE) {
    w_ratio *= lb_new.z * exp(-lb_new.z * k.z) / (2.0 * sinh(lb_new.z * max));
    w_ratio /= lb_old.z * exp(lb_old.z * k.z) / (2.0 * sinh(lb_old.z * max));
  }

  return w_ratio;
}

inline double MultiParticle::GetCoeff()
{
  // calculate (w_new->old/w_old->new) and return it.
  double w_ratio = 1.0;
  double lBeta = lambda * BETA;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);

#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(lBeta, r_max, t_max) reduction(*:w_ratio)
#endif
  for(int m = 0; m < moleculeIndex.size(); m++) {
    uint molNumber = moleculeIndex[m];
    if(moveType == mp::MPROTATE) {
      // rotate: lbt_old, lbt_new are lambda * BETA * torque
      XYZ lbt_old = molTorqueRef.Get(molNumber) * lBeta;
      XYZ lbt_new = molTorqueNew.Get(molNumber) * lBeta;
      w_ratio *= CalculateWRatio(lbt_new, lbt_old, r_k.Get(molNumber), r_max);
    } else {
      // displace: lbf_old, lbf_new are lambda * BETA * force
      XYZ lbf_old = (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) *
                    lBeta;
      XYZ lbf_new = (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) *
                    lBeta;
      w_ratio *= CalculateWRatio(lbf_new, lbf_old, t_k.Get(molNumber), t_max);
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
  double MPCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = MPCoeff * uBoltz;
  double pr = prng();
  bool result = (rejectState == mv::fail_state::NO_FAIL) && pr < accept;
  if(result) {
    sysPotRef = sysPotNew;
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    swap(molForceRef, molForceNew);
    swap(atomForceRef, atomForceNew);
    swap(molForceRecRef, molForceRecNew);
    swap(atomForceRecRef, atomForceRecNew);
    swap(molTorqueRef, molTorqueNew);
    //update reciprocal value
    calcEwald->UpdateRecip(bPick);
  } else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookup);
    calcEwald->exgMolCache();
  }

  moveSetRef.UpdateMoveSettingMultiParticle(bPick, result, moveType);

  moveSetRef.Update(mv::MULTIPARTICLE, result, step, bPick);
}

inline XYZ MultiParticle::CalcRandomTransform(XYZ const &lb, double const max, uint molIndex)
{
  XYZ lbmax = lb * max;
  XYZ num;
  if(std::abs(lbmax.x) > MIN_FORCE && std::abs(lbmax.x) < MAX_FORCE) {
    num.x = log(exp(-1.0 * lbmax.x) + 2 * r123wrapper(molIndex * 3 + 0) * sinh(lbmax.x)) / lb.x;
  } else {
    double rr = r123wrapper(molIndex * 3 + 0) * 2.0 - 1.0;
    num.x = max * rr;
  }

  if(std::abs(lbmax.y) > MIN_FORCE && std::abs(lbmax.y) < MAX_FORCE) {
    num.y = log(exp(-1.0 * lbmax.y) + 2 * r123wrapper(molIndex * 3 + 1) * sinh(lbmax.y)) / lb.y;
  } else {
    double rr = r123wrapper(molIndex * 3 + 1) * 2.0 - 1.0;
    num.y = max * rr;
  }

  if(std::abs(lbmax.z) > MIN_FORCE && std::abs(lbmax.z) < MAX_FORCE) {
    num.z = log(exp(-1.0 * lbmax.z) + 2 * r123wrapper(molIndex * 3 + 2) * sinh(lbmax.z)) / lb.z;
  } else {
    double rr = r123wrapper(molIndex * 3 + 2) * 2.0 - 1.0;
    num.z = max * rr;
  }

  if(num.Length() >= boxDimRef.axis.Min(bPick)) {
    std::cout << "Trial Displacement exceed half of the box length in Multiparticle move.\n";
    std::cout << "Trial transform: " << num;
    exit(EXIT_FAILURE);
  } else if (!std::isfinite(num.Length())) {
    std::cout << "Trial Displacement is not a finite number in Multiparticle move.\n";
    std::cout << "Trial transform: " << num;
    exit(EXIT_FAILURE);
  }

  // We can possibly bound them
  return num;
}

inline void MultiParticle::CalculateTrialDistRot()
{
  uint m, molIndex;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);
  XYZ lbf; // lambda * BETA * force * maxTranslate
  XYZ lbt; // lambda * BETA * torque * maxRotation

  for(m = 0; m < moleculeIndex.size(); m++) {
    molIndex = moleculeIndex[m];

    if(moveType == mp::MPROTATE) { // rotate
      lbt = molTorqueRef.Get(molIndex) * lambda * BETA;
      r_k.Set(molIndex, CalcRandomTransform(lbt, r_max, molIndex));
    } else { // displace
      lbf = (molForceRef.Get(molIndex) + molForceRecRef.Get(molIndex)) *
            lambda * BETA;
      t_k.Set(molIndex, CalcRandomTransform(lbf, t_max, molIndex));
    }
  }
}

inline void MultiParticle::RotateForceBiased(uint molIndex)
{
  XYZ rot = r_k.Get(molIndex);
  double rotLen = rot.Length();
  RotationMatrix matrix;

  XYZ axis = rot * (1.0 / rotLen);
  TransformMatrix cross = TransformMatrix::CrossProduct(axis);
  TransformMatrix tensor = TransformMatrix::TensorProduct(axis);
  matrix = RotationMatrix::FromAxisAngle(rotLen, cross, tensor);

  XYZ center = newCOMs.Get(molIndex);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);

  // Do Rotation
  for(uint p = 0; p < len; p++) {
    temp.Add(p, -center);
    XYZ newPosition = matrix.Apply(temp[p]);
    temp.Set(p, newPosition);
    temp.Add(p, center);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticle::TranslateForceBiased(uint molIndex)
{
  XYZ shift = t_k.Get(molIndex);
  if(shift > boxDimRef.GetHalfAxis(bPick)) {
    std::cout << "Trial Displacement exceed half of the box length in Multiparticle move.\n";
    std::cout << "Trial transform: " << shift;
    exit(EXIT_FAILURE);
  }

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
