/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MULTIPARTICLE_H
#define MULTIPARTICLE_H

#include <cmath>
#include <fstream>

#include "MoveBase.h"
#include "Random123Wrapper.h"
#include "StaticVals.h"
#include "System.h"
#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

class MultiParticle : public MoveBase {
public:
  MultiParticle(System &sys, StaticVals const &statV);
  ~MultiParticle() {}

  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const ulong step);
  virtual void PrintAcceptKind();

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
  bool allTranslate;
  bool multiParticleLiquid, multiParticleGas;
  std::vector<uint> moleculeIndex;
  std::vector<int> inForceRange;
  const MoleculeLookup &molLookup;
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars;
#endif
  Random123Wrapper &r123wrapper;
  const Molecules &mols;

  double GetCoeff();
  uint ChooseBox();
  void CalculateTrialDistRot();
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
  void RotateRandom(uint molIndex);
  void TranslateRandom(uint molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(bool &forceInRange, XYZ const &lb, double const max,
                          uint molIndex);
  double CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old, XYZ const &k,
                         double max);
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV)
    : MoveBase(sys, statV),

      newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
      newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol),
      molLookup(sys.molLookup), r123wrapper(sys.r123wrapper), mols(statV.mol) {
  molTorqueRef.Init(sys.com.Count());
  molTorqueNew.Init(sys.com.Count());
  atomForceRecNew.Init(sys.coordinates.Count());
  molForceRecNew.Init(sys.com.Count());

  t_k.Init(sys.com.Count());
  r_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());
  inForceRange.resize(sys.com.Count(), false);

  // set default value for r_max, t_max, and lambda
  // the value of lambda is based on the paper
  lambda = 0.5;
  for (uint b = 0; b < BOX_TOTAL; b++) {
    initMol[b] = false;
  }

  // Check to see if we have only monoatomic molecules or not
  allTranslate = false;
  uint numAtomsPerKind = 0;
  for (uint k = 0; k < molLookup.GetNumKind(); ++k) {
    numAtomsPerKind += molRef.NumAtoms(k);
  }
  // If we have only one atom in each kind, it means all molecules
  // in the system are monoatomic
  allTranslate = (numAtomsPerKind == molLookup.GetNumKind());
  multiParticleLiquid = sys.statV.multiParticleLiquid;
  multiParticleGas = sys.statV.multiParticleGas;
#ifdef GOMC_CUDA
  cudaVars = sys.statV.forcefield.particles->getCUDAVars();
#endif
}

inline void MultiParticle::PrintAcceptKind() {
  printf("%-37s", "% Accepted MultiParticle ");
  for (uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::MULTIPARTICLE));
  }
  std::cout << std::endl;
}

inline void MultiParticle::SetMolInBox(uint box) {
  // NEED to check if atom is not fixed!
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  moleculeIndex.clear();
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
  MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
  while (thisMol != end) {
    // Make sure this molecule is not fixed in its position
    if (!molLookup.IsFix(*thisMol)) {
      moleculeIndex.push_back(*thisMol);
    }
    thisMol++;
  }
#else
  if (!initMol[box]) {
    moleculeIndex.clear();
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    while (thisMol != end) {
      // Make sure this molecule is not fixed in its position
      if (!molLookup.IsFix(*thisMol)) {
        moleculeIndex.push_back(*thisMol);
      }
      thisMol++;
    }
  }
#endif
  initMol[box] = true;
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MULTIPARTICLE);
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#elif ENSEMBLE == GEMC
  if (multiParticleLiquid != multiParticleGas) // Only pick one of the two boxes
    bPick = ChooseBox();
  else
    prng.PickBox(bPick, subDraw, movPerc);
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif

  // In each step, we perform either:
  // 1- All displacement move.
  // 2- All rotation move.
  if (allTranslate) {
    moveType = mp::MPDISPLACE;
  } else {
    moveType = prng.randIntExc(mp::MPTOTALTYPES);
  }
#ifndef NDEBUG
  if (moveType == mp::MPDISPLACE)
    std::cout << "   MultiParticle Displacement" << std::endl;
  else if (moveType == mp::MPROTATE)
    std::cout << "   MultiParticle Rotation" << std::endl;
  else
    std::cout
        << "   MultiParticle move type not recognized! Update MultiParticle.h"
        << std::endl;
#endif

  SetMolInBox(bPick);
  // reset all inForceRange vector to false.
  std::fill(inForceRange.begin(), inForceRange.end(), false);
  if (moleculeIndex.size() == 0) {
    std::cout << "Warning: MultiParticle move has no particles to move. "
                 "Skipping...\n";
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    return state;
  }

  // We don't use forces for non-MP moves, so we need to calculate them for the
  // current system if any other moves, besides other MP moves, have been
  // accepted. Or, if this is the first MP move, which is handled with the same
  // flag.
  if (moveSetRef.GetSingleMoveAccepted(bPick)) {
    GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
    // Copy ref reciprocal terms to new for calculation with old positions
    calcEwald->CopyRecip(bPick);

    // Calculate long range electrostatic force for old positions
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                  bPick);

    // Calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);

    // Calculate Torque for old positions
    calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                              atomForceRef, atomForceRecRef, molTorqueRef,
                              bPick);

    sysPotRef.Total();
    GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  atomForceRef.CopyRange(atomForceNew, 0, 0, atomForceNew.Count());
  molForceRef.CopyRange(molForceNew, 0, 0, molForceNew.Count());
  atomForceRecRef.CopyRange(atomForceRecNew, 0, 0, atomForceRecNew.Count());
  molForceRecRef.CopyRange(molForceRecNew, 0, 0, molForceRecNew.Count());
  molTorqueRef.CopyRange(molTorqueNew, 0, 0, molTorqueNew.Count());
#endif

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MULTIPARTICLE);
  return state;
}

// To relax the system in NE_MTMC move
inline uint MultiParticle::PrepNEMTMC(const uint box, const uint midx,
                                      const uint kidx) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MULTIPARTICLE);
  uint state = mv::fail_state::NO_FAIL;
  bPick = box;
  // In each step, we perform either:
  // 1- All displacement move.
  // 2- All rotation move.
  if (allTranslate) {
    moveType = mp::MPDISPLACE;
  } else {
    moveType = prng.randIntExc(mp::MPTOTALTYPES);
  }

  SetMolInBox(bPick);
  if (moleculeIndex.size() == 0) {
    std::cout << "Warning: MultiParticle move has no particles to move. "
                 "Skipping...\n";
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    return state;
  }

  // We don't use forces for non-MP moves, so we need to calculate them for the
  // current system if any other moves, besides other MP moves, have been
  // accepted. Or, if this is the first MP move, which is handled with the same
  // flag.
  if (moveSetRef.GetSingleMoveAccepted(bPick)) {
    // Copy ref reciprocal terms to new for calculation with old positions
    calcEwald->CopyRecip(bPick);

    // Calculate long range electrostatic force for old positions
    calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                  bPick);

    // Calculate short range energy and force for old positions
    calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                       boxDimRef, bPick);

    // Calculate Torque for old positions
    calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                              atomForceRef, atomForceRecRef, molTorqueRef,
                              bPick);

    sysPotRef.Total();
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  atomForceRef.CopyRange(atomForceNew, 0, 0, atomForceNew.Count());
  molForceRef.CopyRange(molForceNew, 0, 0, molForceNew.Count());
  atomForceRecRef.CopyRange(atomForceRecNew, 0, 0, atomForceRecNew.Count());
  molForceRecRef.CopyRange(molForceRecNew, 0, 0, molForceRecNew.Count());
  molTorqueRef.CopyRange(molTorqueNew, 0, 0, molTorqueNew.Count());
#endif
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MULTIPARTICLE);
  return state;
}

inline uint MultiParticle::ChooseBox() {
  double minDens = 1.0e20;
  double maxDens = 0.0;
  uint minB = 0, maxB = 0;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    double density = 0.0;
    for (uint k = 0; k < molLookup.GetNumKind(); ++k) {
      density += molLookup.NumKindInBox(k, b) * boxDimRef.volInv[b] *
                 molRef.kinds[k].molMass;
    }
    if (density > maxDens) {
      maxDens = density;
      maxB = b;
    }
    if (density < minDens) {
      minDens = density;
      minB = b;
    }
  }
  if (multiParticleLiquid)
    return maxB;
  else
    return minB;
}

inline uint MultiParticle::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MULTIPARTICLE);
  // Based on the reference force decide whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;
#ifdef GOMC_CUDA
  // One CoM per molecule, so use this to set the vector size.
  std::vector<int8_t> isMoleculeInvolved(newCOMs.Count(), 0);
  // find the IDs of molecules in moleculeIndex
  for (int m = 0; m < (int)moleculeIndex.size(); m++) {
    int mol = moleculeIndex[m];
    isMoleculeInvolved[mol] = 1;
  }

  // This kernel will calculate translation/rotation amount + shifting/rotating
  if (moveType == mp::MPROTATE) {
    double r_max = moveSetRef.GetRMAX(bPick);
    CallRotateParticlesGPU(
        cudaVars, isMoleculeInvolved, bPick, r_max, molTorqueRef.x,
        molTorqueRef.y, molTorqueRef.z, inForceRange, r123wrapper.GetStep(),
        r123wrapper.GetKeyValue(), r123wrapper.GetSeedValue(),
        newMolsPos.Count(), newCOMs.Count(), boxDimRef.GetAxis(bPick).x,
        boxDimRef.GetAxis(bPick).y, boxDimRef.GetAxis(bPick).z, newMolsPos,
        newCOMs, lambda * BETA, r_k);
  } else {
    double t_max = moveSetRef.GetTMAX(bPick);
    CallTranslateParticlesGPU(
        cudaVars, isMoleculeInvolved, bPick, t_max, molForceRef.x,
        molForceRef.y, molForceRef.z, inForceRange, r123wrapper.GetStep(),
        r123wrapper.GetKeyValue(), r123wrapper.GetSeedValue(),
        newMolsPos.Count(), newCOMs.Count(), boxDimRef.GetAxis(bPick).x,
        boxDimRef.GetAxis(bPick).y, boxDimRef.GetAxis(bPick).z, newMolsPos,
        newCOMs, lambda * BETA, t_k, molForceRecRef);
  }
#else
  // Calculate trial translate and rotate
  // move particles according to force and torque and store them in the new pos
  CalculateTrialDistRot();
#endif

  // Do error checking and skip if there is an invalid transform amount.
  for (int m = 0; m < (int)moleculeIndex.size(); m++) {
    uint molIndex = moleculeIndex[m];
    XYZ num;
    if (moveType == mp::MPROTATE) // rotate
      num = r_k.Get(molIndex);
    else { // displace
      num = t_k.Get(molIndex);
      // check for PBC error and bad initial configuration
      if (num > boxDimRef.GetHalfAxis(bPick)) {
        std::cout << "Warning: Trial Displacement exceeds half the box length "
                     "in Multiparticle move!"
                  << std::endl;
        std::cout << "         Trial transformation vector: " << num
                  << std::endl;
        std::cout << "         Box Dimensions: " << boxDimRef.GetAxis(bPick)
                  << std::endl
                  << std::endl;
        std::cout << "This might be due to a bad initial configuration, where "
                     "atoms of the molecules"
                  << std::endl
                  << "are too close to each other or overlap. Please "
                     "equilibrate your system using"
                  << std::endl
                  << "rigid body translation or rotation MC moves before using "
                     "the Multiparticle"
                  << std::endl
                  << "move." << std::endl
                  << std::endl;
        state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
        break;
      }
    }
    if (!std::isfinite(num.Length())) {
      std::cout << "Trial Displacement is not a finite number in MultiParticle";
      std::cout << " move.\nTrial transform: " << num;
      state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
      break;
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MULTIPARTICLE);
  return state;
}

inline void MultiParticle::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  cellList.GridAll(boxDimRef, newMolsPos, molLookup);

  // back up cached Fourier term
  calcEwald->backupMolCache();

  // setup reciprocal vectors for new positions
  calcEwald->BoxReciprocalSums(bPick, newMolsPos);

  sysPotNew = sysPotRef;
  // calculate short range energy and force
  sysPotNew = calcEnRef.BoxForce(sysPotNew, newMolsPos, atomForceNew,
                                 molForceNew, boxDimRef, bPick);
  // calculate long range electrostatic energy for new positions
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick, false);
  // Calculate long range electrostatic force for new positions
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecNew, molForceRecNew,
                                bPick);

  // Calculate Torque for new positions
  calcEnRef.CalculateTorque(moleculeIndex, newMolsPos, newCOMs, atomForceNew,
                            atomForceRecNew, molTorqueNew, bPick);

  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
}

inline double MultiParticle::CalculateWRatio(XYZ const &lb_new,
                                             XYZ const &lb_old, XYZ const &k,
                                             double max) {
  double w_ratio = 1.0;

  w_ratio *= lb_new.x * exp(-lb_new.x * k.x) / (2.0 * sinh(lb_new.x * max));
  w_ratio /= lb_old.x * exp(lb_old.x * k.x) / (2.0 * sinh(lb_old.x * max));

  w_ratio *= lb_new.y * exp(-lb_new.y * k.y) / (2.0 * sinh(lb_new.y * max));
  w_ratio /= lb_old.y * exp(lb_old.y * k.y) / (2.0 * sinh(lb_old.y * max));

  w_ratio *= lb_new.z * exp(-lb_new.z * k.z) / (2.0 * sinh(lb_new.z * max));
  w_ratio /= lb_old.z * exp(lb_old.z * k.z) / (2.0 * sinh(lb_old.z * max));

  return w_ratio;
}

inline double MultiParticle::GetCoeff() {
  // calculate (w_new->old/w_old->new) and return it.
  double w_ratio = 1.0;
  double lBeta = lambda * BETA;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);

  if (moveType == mp::MPROTATE) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(lBeta, r_max, t_max) reduction(*:w_ratio)
#endif
    for (int m = 0; m < (int)moleculeIndex.size(); m++) {
      uint molNumber = moleculeIndex[m];
      // If force or torque was not in the range, no need to calculate weight
      // ratio it's simply 1.0
      if (inForceRange[molNumber]) {
        // rotate: lbt_old, lbt_new are lambda * BETA * torque
        XYZ lbt_old = molTorqueRef.Get(molNumber) * lBeta;
        XYZ lbt_new = molTorqueNew.Get(molNumber) * lBeta;
        w_ratio *= CalculateWRatio(lbt_new, lbt_old, r_k.Get(molNumber), r_max);
      }
    }
  } else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(lBeta, r_max, t_max) reduction(*:w_ratio)
#endif
    for (int m = 0; m < (int)moleculeIndex.size(); m++) {
      uint molNumber = moleculeIndex[m];
      // If force or torque was not in the range, no need to calculate weight
      // ratio it's simply 1.0
      if (inForceRange[molNumber]) {
        // displace: lbf_old, lbf_new are lambda * BETA * force
        XYZ lbf_old =
            (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) *
            lBeta;
        XYZ lbf_new =
            (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) *
            lBeta;
        w_ratio *= CalculateWRatio(lbf_new, lbf_old, t_k.Get(molNumber), t_max);
      }
    }
  }

  // In case where force or torque is a large negative number (ex. -800) the
  // exp value becomes inf. In these situations we return 0.0 to reject the move
  if (!std::isfinite(w_ratio)) {
    w_ratio = 0.0;
  }

  return w_ratio;
}

inline void MultiParticle::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_MULTIPARTICLE);
  // Here we compare the values of reference and trial and decide whether to
  // accept or reject the move
  double MPCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = MPCoeff * uBoltz;
  double pr = r123wrapper.GetRandomNumber(newCOMs.Count());
  bool result = (rejectState == mv::fail_state::NO_FAIL) && pr < accept;
  if (result) {
    sysPotRef = sysPotNew;
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    swap(molForceRef, molForceNew);
    swap(atomForceRef, atomForceNew);
    swap(molForceRecRef, molForceRecNew);
    swap(atomForceRecRef, atomForceRecNew);
    swap(molTorqueRef, molTorqueNew);
    // Update reciprocal value
    calcEwald->UpdateRecip(bPick);
    // Update the velocity in box
    velocity.UpdateBoxVelocity(bPick);
  } else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookup);
    calcEwald->exgMolCache();
  }

  moveSetRef.UpdateMoveSettingMultiParticle(bPick, result, moveType);
  moveSetRef.Update(mv::MULTIPARTICLE, result, bPick);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_MULTIPARTICLE);
}

inline XYZ MultiParticle::CalcRandomTransform(bool &forceInRange, XYZ const &lb,
                                              double const max, uint molIndex) {
  XYZ lbmax = lb * max;
  // XYZ default constructor initializes to (0.0, 0.0, 0.0)
  XYZ val;
  forceInRange =
      std::abs(lbmax.x) > MIN_FORCE && std::abs(lbmax.x) < MAX_FORCE &&
      std::abs(lbmax.y) > MIN_FORCE && std::abs(lbmax.y) < MAX_FORCE &&
      std::abs(lbmax.z) > MIN_FORCE && std::abs(lbmax.z) < MAX_FORCE;

  if (forceInRange) {
    XYZ randnums = r123wrapper.GetRandomCoords(molIndex);
    val.x = log(exp(-1.0 * lbmax.x) + 2.0 * randnums.x * sinh(lbmax.x)) / lb.x;
    val.y = log(exp(-1.0 * lbmax.y) + 2.0 * randnums.y * sinh(lbmax.y)) / lb.y;
    val.z = log(exp(-1.0 * lbmax.z) + 2.0 * randnums.z * sinh(lbmax.z)) / lb.z;
  }

  return val;
}

inline void MultiParticle::CalculateTrialDistRot() {
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);

  if (moveType == mp::MPROTATE) { // rotate
    double *x = r_k.x;
    double *y = r_k.y;
    double *z = r_k.z;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(r_max, x, y, z)
#endif
    for (uint m = 0; m < moleculeIndex.size(); m++) {
      uint molIndex = moleculeIndex[m];
      XYZ lbt = molTorqueRef.Get(molIndex) * lambda * BETA;
      bool forceInRange;
      XYZ val = CalcRandomTransform(forceInRange, lbt, r_max, molIndex);
      x[molIndex] = val.x;
      y[molIndex] = val.y;
      z[molIndex] = val.z;
      inForceRange[molIndex] = forceInRange;
      if (forceInRange) {
        RotateForceBiased(molIndex);
      } else {
        RotateRandom(molIndex);
      }
    }
  } else { // displace
    double *x = t_k.x;
    double *y = t_k.y;
    double *z = t_k.z;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(t_max, x, y, z)
#endif
    for (uint m = 0; m < moleculeIndex.size(); m++) {
      uint molIndex = moleculeIndex[m];
      XYZ lbf = (molForceRef.Get(molIndex) + molForceRecRef.Get(molIndex)) *
                lambda * BETA;
      bool forceInRange;
      XYZ val = CalcRandomTransform(forceInRange, lbf, t_max, molIndex);
      x[molIndex] = val.x;
      y[molIndex] = val.y;
      z[molIndex] = val.z;
      inForceRange[molIndex] = forceInRange;
      if (forceInRange) {
        TranslateForceBiased(molIndex);
      } else {
        TranslateRandom(molIndex);
      }
    }
  }
}

inline void MultiParticle::RotateForceBiased(uint molIndex) {
  XYZ rot = r_k.Get(molIndex);
  double rotLen = rot.Length();
  RotationMatrix matrix;

  XYZ axis = rot * (1.0 / rotLen);
  TransformMatrix cross = TransformMatrix::CrossProduct(axis);
  TransformMatrix tensor = TransformMatrix::TensorProduct(axis);
  matrix = RotationMatrix::FromAxisAngle(rotLen, cross, tensor);

  XYZ center = comCurrRef.Get(molIndex);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);

  // Do Rotation
  for (uint p = 0; p < len; p++) {
    temp.Add(p, -center);
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, center);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticle::TranslateForceBiased(uint molIndex) {
  XYZ shift = t_k.Get(molIndex);
  XYZ newcom = comCurrRef.Get(molIndex);
  uint stop, start, len;
  molRef.GetRange(start, stop, len, molIndex);
  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  // Shift the coordinate and COM
  temp.AddAll(shift);
  newcom += shift;
  // rewrapping
  boxDimRef.WrapPBC(temp, bPick);
  newcom = boxDimRef.WrapPBC(newcom, bPick);
  // set the new coordinate
  temp.CopyRange(newMolsPos, 0, start, len);
  newCOMs.Set(molIndex, newcom);
}

inline void MultiParticle::RotateRandom(uint molIndex) {
  double r_max = moveSetRef.GetRMAX(bPick);
  double symRand = r123wrapper.GetSymRandom(molIndex, r_max);
  XYZ sphereCoords = r123wrapper.GetRandomCoordsOnSphere(molIndex);
  RotationMatrix matrix = RotationMatrix::FromAxisAngle(symRand, sphereCoords);

  XYZ center = comCurrRef.Get(molIndex);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);

  // Do Rotation
  for (uint p = 0; p < len; p++) {
    temp.Add(p, -center);
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, center);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticle::TranslateRandom(uint molIndex) {
  double t_max = moveSetRef.GetTMAX(bPick);
  XYZ shift = r123wrapper.GetSymRandomCoords(molIndex, t_max);
  XYZ newcom = comCurrRef.Get(molIndex);
  uint stop, start, len;

  molRef.GetRange(start, stop, len, molIndex);
  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  // Shift the coordinate and COM
  temp.AddAll(shift);
  newcom += shift;
  // rewrapping
  boxDimRef.WrapPBC(temp, bPick);
  newcom = boxDimRef.WrapPBC(newcom, bPick);
  // set the new coordinate
  temp.CopyRange(newMolsPos, 0, start, len);
  newCOMs.Set(molIndex, newcom);
}

#endif
