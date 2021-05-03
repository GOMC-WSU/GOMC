/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MULTIPARTICLEBROWNIANMOTION_H
#define MULTIPARTICLEBROWNIANMOTION_H

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"
#include <cmath>
#include "Random123Wrapper.h"
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "CUDAMemoryManager.cuh"
#include "TransformParticlesCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif

class MultiParticleBrownian : public MoveBase
{
public:
  MultiParticleBrownian(System &sys, StaticVals const& statV);
  ~MultiParticleBrownian() {
#ifdef GOMC_CUDA
    cudaVars = NULL;
    cudaFreeHost(kill);
    kill = NULL;
#endif
  }

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();
  //used in CFCMC for initialization
  void PrepCFCMC(const uint box);

private:
  uint bPick;
  bool initMol;
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
  Random123Wrapper &r123wrapper;
  bool allTranslate;
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars;
  bool isOrthogonal;
  int *kill; // kill the simulation if we started with bad configuration
#endif

  double GetCoeff();
  void CalculateTrialDistRot();
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(XYZ const &lb, double const max, uint molIndex);
  double CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old, XYZ const &k,
                         double max4);
};

inline MultiParticleBrownian::MultiParticleBrownian(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV),
  newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol),
  molLookup(sys.molLookup), r123wrapper(sys.r123wrapper)
{
  molTorqueNew.Init(sys.com.Count());
  molTorqueRef.Init(sys.com.Count());
  atomForceRecNew.Init(sys.coordinates.Count());
  molForceRecNew.Init(sys.com.Count());

  t_k.Init(sys.com.Count());
  r_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());

  initMol = false;
  
  // Check to see if we have only monoatomic molecule or not
  allTranslate = false;
  int numAtomsPerKind = 0;
  for (int k = 0; k < molLookup.GetNumKind(); ++k) {
    numAtomsPerKind += molRef.NumAtoms(k);
  }
  // If we have only one atom in each kind, it means all molecule
  // in the system is monoatomic
  allTranslate = (numAtomsPerKind == molLookup.GetNumKind());

#ifdef GOMC_CUDA
  cudaVars = sys.statV.forcefield.particles->getCUDAVars();
  isOrthogonal = statV.isOrthogonal;
  cudaMallocHost((void**) &kill, sizeof(int));
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

inline void MultiParticleBrownian::PrintAcceptKind()
{
  printf("%-37s", "% Accepted MultiParticle-Brownian ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::MULTIPARTICLE_BM));
  }
  std::cout << std::endl;
}


inline void MultiParticleBrownian::SetMolInBox(uint box)
{
  // NEED to check if atom is not fixed!
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  // need to be initialized for every move since number of atom is changing
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
  // box would be always 0 in NVT or NPT ensemble, initialize it once
  if(!initMol) {
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
    initMol = true;
  }
#endif
}

inline uint MultiParticleBrownian::Prep(const double subDraw, const double movPerc)
{
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif

  // In each step, we perform either:
  // 1- All displacement move.
  // 2- All rotation move.
  if(allTranslate) {
    moveType = mp::MPDISPLACE;
  } else {
    moveType = prng.randIntExc(mp::MPTOTALTYPES);
  }

  SetMolInBox(bPick);
  if (moleculeIndex.size() == 0) {
    std::cout << "Warning: MultiParticleBrownian move can't move any molecules, Skipping...\n";
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    return state;
  }

  //We don't use forces for non-MP moves, so we need to calculate them for the
  //current system if any other moves, besides other MP moves, have been accepted.
  //Or, if this is the first MP move, which is handled with the same flag.
  if(moveSetRef.GetSingleMoveAccepted(bPick)) {
    GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);
    //Copy ref reciprocal terms to new for calculation with old positions
    calcEwald->CopyRecip(bPick);

    //Calculate force for long range electrostatic using old position
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef, bPick);

    //calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);

    if(moveType == mp::MPROTATE) {
      //Calculate Torque for old positions
      calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                                atomForceRef, atomForceRecRef, molTorqueRef, bPick);
    }

    sysPotRef.Total();
    GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
  return state;
}

inline void MultiParticleBrownian::PrepCFCMC(const uint box)
{
  bPick = box;
  // In each step, we perform either:
  // 1- All displacement move.
  // 2- All rotation move.
  if(allTranslate) {
    moveType = mp::MPDISPLACE;
  } else {
    moveType = prng.randIntExc(mp::MPTOTALTYPES);
  }

  SetMolInBox(bPick);
  //We don't use forces for non-MP moves, so we need to calculate them for the
  //current system if any other moves, besides other MP moves, have been accepted.
  //Or, if this is the first MP move, which is handled with the same flag.
  if(moveSetRef.GetSingleMoveAccepted(bPick)) {
    //Copy ref reciprocal terms to new for calculation with old positions
    calcEwald->CopyRecip(bPick);

    //Calculate long range electrostatic force for old positions
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef, bPick);

    //Calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);

    if(moveType == mp::MPROTATE) {
      //Calculate Torque for old positions
      calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                                atomForceRef, atomForceRecRef, molTorqueRef, bPick);
    }

    sysPotRef.Total();
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
}

inline uint MultiParticleBrownian::Transform()
{
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MULTIPARTICLE_BM);
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;

#ifdef GOMC_CUDA
  kill[0] = 0;
  // This kernel will calculate translation/rotation amount + shifting/rotating
  if(moveType == mp::MPROTATE) {
    double r_max = moveSetRef.GetRMAX(bPick);
    BrownianMotionRotateParticlesGPU(
      cudaVars,
      moleculeIndex,
      molTorqueRef,
      newMolsPos, 
      newCOMs, 
      r_k, 
      boxDimRef.GetAxis(bPick),
      BETA, 
      r_max,
      r123wrapper.GetStep(), 
      r123wrapper.GetSeedValue(),
      bPick,
      isOrthogonal,
      kill);
  } else {
    double t_max = moveSetRef.GetTMAX(bPick);
    BrownianMotionTranslateParticlesGPU(
      cudaVars,
      moleculeIndex,
      molForceRef,
      molForceRecRef,
      newMolsPos, 
      newCOMs, 
      t_k, 
      boxDimRef.GetAxis(bPick),
      BETA, 
      t_max,
      r123wrapper.GetStep(), 
      r123wrapper.GetSeedValue(),
      bPick,
      isOrthogonal,
      kill);
  }
  // kill the simulation if we had bad initial configuration
  if(kill[0]) {
    std::cout << "Error: Transformation of " << kill[0] << "in Multiparticle Brownian Motion move failed!\n"
              << "       Either trial rotate/translate is not finite number or their translation \n" 
              << "       exceeded half of the box length! \n\n";
    std::cout << "This might be due to bad initial configuration, where atom of the molecules \n" 
              << "are too close to each other or overlaps. Please equilibrate your system using \n"
              << "rigid body translation or rotation MC moves, before using the Multiparticle \n"
              << "Brownian Motion move. \n\n";
    exit(EXIT_FAILURE);
} 
#else
  // Calculate trial translate and rotate
  // move particles according to force and torque and store them in the new pos
  CalculateTrialDistRot();
#endif
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MULTIPARTICLE_BM);
  return state;
}

inline void MultiParticleBrownian::CalcEn()
{
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  //cellList.GridAll(boxDimRef, newMolsPos, molLookup);
  cellList.GridBox(boxDimRef, newMolsPos, molLookup, bPick);

  //back up cached fourier term
  calcEwald->backupMolCache();

  //setup reciprocal vectors for new positions
  calcEwald->BoxReciprocalSums(bPick, newMolsPos);

  sysPotNew = sysPotRef;
  //calculate short range energy and force
  sysPotNew = calcEnRef.BoxForce(sysPotNew, newMolsPos, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  //calculate long range of new electrostatic energy
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick, false);
  //Calculate long range of new electrostatic force
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecNew, molForceRecNew,
                                bPick);

  if(moveType == mp::MPROTATE) {
    //Calculate Torque for new positions
    calcEnRef.CalculateTorque(moleculeIndex, newMolsPos, newCOMs, atomForceNew,
                              atomForceRecNew, molTorqueNew, bPick);
  }
  sysPotNew.Total();
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
}

inline double MultiParticleBrownian::CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old,
    XYZ const &k, double max4)
{
  double w_ratio = 0.0;
  XYZ old_var = lb_old - k;
  XYZ new_var = lb_new + k;

  //Note: we could factor max4 and multiply at the end, but
  //      for the move, where we translate and rotate all molecules,
  //      this method would not work. Hence, I did not factor it.
  // its actually is w_ratio += -1.0* but we simplify it
  w_ratio -= (new_var.LengthSq() / max4);
  // its actually is w_ratio -= -1.0* but we simplify it
  w_ratio += (old_var.LengthSq() / max4);

  return w_ratio;
}

inline double MultiParticleBrownian::GetCoeff()
{
  // calculate (w_new->old/w_old->new) and return it.
  double w_ratio = 0.0;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);
  double r_max4 = 4.0 * r_max;
  double t_max4 = 4.0 * t_max;

  if(moveType == mp::MPROTATE) {// rotate, 
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(moleculeIndex, r_max, t_max, r_max4, t_max4) reduction(+:w_ratio)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      uint molNumber = moleculeIndex[m];
      // rotate, bt_ = BETA * force * maxForce
      XYZ bt_old = molTorqueRef.Get(molNumber) * BETA * r_max;
      XYZ bt_new = molTorqueNew.Get(molNumber) * BETA * r_max;
      w_ratio += CalculateWRatio(bt_new, bt_old, r_k.Get(molNumber), r_max4);
    } 
  } else {// displace
#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(moleculeIndex, r_max, t_max, r_max4, t_max4) reduction(+:w_ratio)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      uint molNumber = moleculeIndex[m];
      // bf_ = BETA * torque * maxTorque
      XYZ bf_old = (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) *
                BETA * t_max;
      XYZ bf_new = (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) *
                BETA * t_max;
      w_ratio += CalculateWRatio(bf_new, bf_old, t_k.Get(molNumber), t_max4);
    }
  }

  return w_ratio;
}

inline void MultiParticleBrownian::Accept(const uint rejectState, const uint step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_MULTIPARTICLE_BM);
  // Here we compare the values of reference and trial and decide whether to
  // accept or reject the move
  double MPCoeff = GetCoeff();
  double accept = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()) + MPCoeff);
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
    // Update the velocity in box
    velocity.UpdateBoxVelocity(bPick);
  } else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookup);
    calcEwald->exgMolCache();
  }

  moveSetRef.UpdateMoveSettingMultiParticle(bPick, result, moveType);
  moveSetRef.Update(mv::MULTIPARTICLE_BM, result, bPick);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_MULTIPARTICLE_BM);
}

inline XYZ MultiParticleBrownian::CalcRandomTransform(XYZ const &lb, double const max, uint molIndex)
{
  XYZ lbmax = lb * max;
  XYZ num;
  //variance is 2A according to the paper, so stdDev is sqrt(variance)
  double stdDev = sqrt(2.0 * max);

  num.x = lbmax.x + r123wrapper.GetGaussianNumber(molIndex * 3 + 0, 0.0, stdDev);
  num.y = lbmax.y + r123wrapper.GetGaussianNumber(molIndex * 3 + 1, 0.0, stdDev);
  num.z = lbmax.z + r123wrapper.GetGaussianNumber(molIndex * 3 + 2, 0.0, stdDev);


  if (!std::isfinite(num.Length())) {
    std::cout << "Error: Trial transform is not a finite number in Brownian Motion Multiparticle move.\n";
    std::cout << "       Trial transform: " << num;
    exit(EXIT_FAILURE);
  }
  // We can possible bound them
  return num;
}

inline void MultiParticleBrownian::CalculateTrialDistRot()
{
  uint molIndex;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);
  if(moveType == mp::MPROTATE) { // rotate
    double *x = r_k.x;
    double *y = r_k.y;
    double *z = r_k.z;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(moleculeIndex, r_max, x, y, z)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      uint molIndex = moleculeIndex[m];
      XYZ bt = molTorqueRef.Get(molIndex) * BETA;
      XYZ val = CalcRandomTransform(bt, r_max, molIndex);
      x[molIndex] = val.x; 
      y[molIndex] = val.y; 
      z[molIndex] = val.z;
      RotateForceBiased(molIndex);
    }
  } else { // displace
    double *x = t_k.x;
    double *y = t_k.y;
    double *z = t_k.z;
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(moleculeIndex, t_max, x, y, z)
#endif
    for(int m = 0; m < moleculeIndex.size(); m++) {
      uint molIndex = moleculeIndex[m];
      XYZ bf = (molForceRef.Get(molIndex) + molForceRecRef.Get(molIndex)) * BETA;
      XYZ val = CalcRandomTransform(bf, t_max, molIndex);
      x[molIndex] = val.x; 
      y[molIndex] = val.y; 
      z[molIndex] = val.z; 
      TranslateForceBiased(molIndex);
    }
  }
}

inline void MultiParticleBrownian::RotateForceBiased(uint molIndex)
{
  XYZ rot = r_k.Get(molIndex);
  double rotLen = rot.Length();
  RotationMatrix matrix;

  matrix = RotationMatrix::FromAxisAngle(rotLen, rot.Normalize());

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
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, center);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticleBrownian::TranslateForceBiased(uint molIndex)
{
  XYZ shift = t_k.Get(molIndex);
  // check for PBC error and bad initial configuration
  if(shift > boxDimRef.GetHalfAxis(bPick)) {
    std::cout << "Error: Trial Displacement exceed half of the box length in Multiparticle \n" 
              << "       Brownian Motion move!\n";
    std::cout << "       Trial transformation vector: " << shift << std::endl;
    std::cout << "       Box Dimension: " << boxDimRef.GetAxis(bPick) << std::endl << std::endl;
    std::cout << "This might be due to bad initial configuration, where atom of the molecules \n" 
              << "are too close to each other or overlaps. Please equilibrate your system using \n"
              << "rigid body translation or rotation MC moves, before using the Multiparticle \n"
              << "Brownian Motion move. \n\n";
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
