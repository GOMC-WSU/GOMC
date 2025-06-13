/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MULTIPARTICLE_H
#define MULTIPARTICLE_H

#include "MoveBase.h"
#include "Random123Wrapper.h"
#include "StaticVals.h"
#include "System.h"
#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif
#include <fstream>

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

class MultiParticle : public MoveBase {
public:
  MultiParticle(System &sys, StaticVals const &statV);
  virtual ~MultiParticle() {}

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
  bool initMol; // Used for ensembles with only one box
  SystemPotential sysPotNew;
  XYZArray molForceNew;
  XYZArray molForceRecNew;
  XYZArray molTorqueRef;
  XYZArray molTorqueNew;
  XYZArray rt_k;
  Coordinates newMolsPos;
  COM newCOMs;
  int moveType;
  bool allTranslate;
  bool multiParticleLiquid, multiParticleGas;
  std::vector<uint> moleculeIndex;
  const MoleculeLookup &molLookup;
  std::vector<int8_t> inForceRange;
  Random123Wrapper &r123wrapper;
  const Molecules &mols;
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars;
#endif

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
  molForceNew.Init(sys.com.Count());
  molForceRecNew.Init(sys.com.Count());
  molTorqueRef.Init(sys.com.Count());
  molTorqueNew.Init(sys.com.Count());

  rt_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());
  inForceRange.resize(sys.com.Count(), false);

  // set default value for r_max, t_max, and lambda
  // the value of lambda is based on the paper
  lambda = 0.5;
  initMol = false;

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
  // box is always 0 in NVT or NPT ensemble, initialize it once
  if (!initMol) {
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
    initMol = true;
  }
#endif
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
    std::cout << "   MultiParticle Displacement\n";
  else if (moveType == mp::MPROTATE)
    std::cout << "   MultiParticle Rotation\n";
  else
    std::cout << "   MultiParticle move type not recognized! Update "
                 "MultiParticle.h\n";
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

  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
  // Copy ref reciprocal terms to new for calculation with old positions
  calcEwald->CopyRecip(bPick);

  // Calculate short range energy and force for old positions
  calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                     boxDimRef, moveType, bPick, false);

  // Calculate long range electrostatic force for old positions
  calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                moveType, bPick);

  if (moveType == mp::MPROTATE) {
    // Calculate Torque for old positions
    calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                              atomForceRef, atomForceRecRef, molTorqueRef,
                              bPick);
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);

  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());

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

  // Copy ref reciprocal terms to new for calculation with old positions
  calcEwald->CopyRecip(bPick);

  // Calculate short range energy and force for old positions
  calcEnRef.BoxForce(sysPotRef, coordCurrRef, atomForceRef, molForceRef,
                     boxDimRef, moveType, bPick, false);

  // Calculate long range electrostatic force for old positions
  calcEwald->BoxForceReciprocal(coordCurrRef, atomForceRecRef, molForceRecRef,
                                moveType, bPick);

  if (moveType == mp::MPROTATE) {
    // Calculate Torque for old positions
    calcEnRef.CalculateTorque(moleculeIndex, coordCurrRef, comCurrRef,
                              atomForceRef, atomForceRecRef, molTorqueRef,
                              bPick);
  }

  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());

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
        newCOMs, lambda * BETA, rt_k);
  } else {
    double t_max = moveSetRef.GetTMAX(bPick);
    CallTranslateParticlesGPU(
        cudaVars, isMoleculeInvolved, bPick, t_max, molForceRef.x,
        molForceRef.y, molForceRef.z, inForceRange, r123wrapper.GetStep(),
        r123wrapper.GetKeyValue(), r123wrapper.GetSeedValue(),
        newMolsPos.Count(), newCOMs.Count(), boxDimRef.GetAxis(bPick).x,
        boxDimRef.GetAxis(bPick).y, boxDimRef.GetAxis(bPick).z, newMolsPos,
        newCOMs, lambda * BETA, rt_k, molForceRecRef);
  }
#else
  // Calculate trial translate and rotate
  // move particles according to force and torque and store them in the new pos
  CalculateTrialDistRot();
#endif

  // Do error checking and skip if there is an invalid transform amount.
  for (int m = 0; m < (int)moleculeIndex.size(); m++) {
    uint molIndex = moleculeIndex[m];
    XYZ num = rt_k.Get(molIndex);
    if (moveType == mp::MPDISPLACE) { // displace
      // check for PBC error and bad initial configuration
      if (num > boxDimRef.GetHalfAxis(bPick)) {
        std::cout << "Warning: Trial Displacement exceeds half the box length "
                     "in Multiparticle move!\n";
        std::cout << "         Trial transformation vector: " << num << "\n";
        std::cout << "         Box Dimensions: " << boxDimRef.GetAxis(bPick);
        std::cout << "\n\n";
        std::cout << "This might be due to a bad initial configuration, where "
                     "atoms of the molecules\n"
                  << "are too close to each other or overlap. Please "
                     "equilibrate your system using\n"
                  << "rigid body translation or rotation MC moves before using "
                     "the Multiparticle move.\n\n";
        state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
        break;
      }
    }
    if (!std::isfinite(num.Length())) {
      std::cout << "Trial Displacement is not a finite number in MultiParticle";
      std::cout << " move.\nTrial transform: " << num << "\n";
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

  // calculate short range energy and force
  // this updates the Lennard-Jones and real electrostatic energy and
  // assigns the result to the new system potential
  sysPotNew = calcEnRef.BoxForce(sysPotRef, newMolsPos, atomForceRef,
                                 molForceNew, boxDimRef, moveType, bPick, true);

  // calculate long range electrostatic energy for new positions
  sysPotNew.boxEnergy[bPick].recip = calcEwald->BoxReciprocal(bPick, false);

  // Calculate long range electrostatic force for new positions
  calcEwald->BoxForceReciprocal(newMolsPos, atomForceRecRef, molForceRecNew,
                                moveType, bPick);

  if (moveType == mp::MPROTATE) {
    // Calculate Torque for new positions
    calcEnRef.CalculateTorque(moleculeIndex, newMolsPos, newCOMs, atomForceRef,
                              atomForceRecRef, molTorqueNew, bPick);
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE);
}

inline double MultiParticle::CalculateWRatio(XYZ const &lb_new,
                                             XYZ const &lb_old, XYZ const &k,
                                             double max) {
  double w_ratio = 1.0;

  w_ratio *=
      lb_new.x * std::exp(-lb_new.x * k.x) / (2.0 * std::sinh(lb_new.x * max));
  w_ratio /=
      lb_old.x * std::exp(lb_old.x * k.x) / (2.0 * std::sinh(lb_old.x * max));

  w_ratio *=
      lb_new.y * std::exp(-lb_new.y * k.y) / (2.0 * std::sinh(lb_new.y * max));
  w_ratio /=
      lb_old.y * std::exp(lb_old.y * k.y) / (2.0 * std::sinh(lb_old.y * max));

  w_ratio *=
      lb_new.z * std::exp(-lb_new.z * k.z) / (2.0 * std::sinh(lb_new.z * max));
  w_ratio /=
      lb_old.z * std::exp(lb_old.z * k.z) / (2.0 * std::sinh(lb_old.z * max));

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
        w_ratio *=
            CalculateWRatio(lbt_new, lbt_old, rt_k.Get(molNumber), r_max);
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
        w_ratio *=
            CalculateWRatio(lbf_new, lbf_old, rt_k.Get(molNumber), t_max);
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
  // Compare the energies of reference and trial positions to decide whether to
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
  XYZ lbmx = lb * max;
  // XYZ default constructor initializes to (0.0, 0.0, 0.0)
  XYZ val;
  forceInRange = std::abs(lbmx.x) > MIN_FORCE && std::abs(lbmx.x) < MAX_FORCE &&
                 std::abs(lbmx.y) > MIN_FORCE && std::abs(lbmx.y) < MAX_FORCE &&
                 std::abs(lbmx.z) > MIN_FORCE && std::abs(lbmx.z) < MAX_FORCE;

  if (forceInRange) {
    XYZ randnums = r123wrapper.GetRandomCoords(molIndex);
    val.x = std::log(std::exp(-lbmx.x) + 2.0 * randnums.x * std::sinh(lbmx.x)) /
            lb.x;
    val.y = std::log(std::exp(-lbmx.y) + 2.0 * randnums.y * std::sinh(lbmx.y)) /
            lb.y;
    val.z = std::log(std::exp(-lbmx.z) + 2.0 * randnums.z * std::sinh(lbmx.z)) /
            lb.z;
  }

  return val;
}

inline void MultiParticle::CalculateTrialDistRot() {
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);

  double *x = rt_k.x;
  double *y = rt_k.y;
  double *z = rt_k.z;

  if (moveType == mp::MPROTATE) { // rotate
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
  XYZ rot = rt_k.Get(molIndex);
  double rotLen = rot.Length();
  RotationMatrix matrix;

  XYZ axis = rot * (1.0 / rotLen);
  TransformMatrix cross = TransformMatrix::CrossProduct(axis);
  TransformMatrix tensor = TransformMatrix::TensorProduct(axis);
  matrix = RotationMatrix::FromAxisAngle(rotLen, cross, tensor);

  XYZ center = comCurrRef.Get(molIndex);
  XYZ CoMoffset = center - matrix.Apply(center);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);

  // Computing the rotation this way has less round off than other methods
  for (uint p = 0; p < len; ++p) { // Rotate each atom in the molecule
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, CoMoffset);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticle::TranslateForceBiased(uint molIndex) {
  XYZ shift = rt_k.Get(molIndex);
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
