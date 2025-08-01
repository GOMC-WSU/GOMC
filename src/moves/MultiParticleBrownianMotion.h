/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MULTIPARTICLEBROWNIANMOTION_H
#define MULTIPARTICLEBROWNIANMOTION_H

#include "MoveBase.h"
#include "Random123Wrapper.h"
#include "StaticVals.h"
#include "System.h"
#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif

class MultiParticleBrownian : public MoveBase {
public:
  MultiParticleBrownian(System &sys, StaticVals const &statV);
  virtual ~MultiParticleBrownian() {}

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
  std::vector<int> moleculeIndex;
  const MoleculeLookup &molLookup;
  Random123Wrapper &r123wrapper;
#ifdef GOMC_CUDA
  VariablesCUDA *cudaVars;
  bool isOrthogonal;
#endif

  double GetCoeff();
  uint ChooseBox();
  void CalculateTrialDistRot();
  void RotateForceBiased(int molIndex);
  void TranslateForceBiased(int molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(XYZ const &lb, double const max, int molIndex);
  double CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old, XYZ const &k,
                         double max4);
};

inline MultiParticleBrownian::MultiParticleBrownian(System &sys,
                                                    StaticVals const &statV)
    : MoveBase(sys, statV),
      newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
      newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol),
      molLookup(sys.molLookup), r123wrapper(sys.r123wrapper) {
  molForceNew.Init(sys.com.Count());
  molForceRecNew.Init(sys.com.Count());
  molTorqueRef.Init(sys.com.Count());
  molTorqueNew.Init(sys.com.Count());

  rt_k.Init(sys.com.Count());
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());

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
  isOrthogonal = statV.isOrthogonal;
#endif
}

inline void MultiParticleBrownian::PrintAcceptKind() {
  printf("%-37s", "% Accepted MultiParticle-Brownian ");
  for (uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::MULTIPARTICLE_BM));
  }
  std::cout << std::endl;
}

inline void MultiParticleBrownian::SetMolInBox(uint box) {
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  // need to be initialized for every move since the number of atoms is changing
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

inline uint MultiParticleBrownian::Prep(const double subDraw,
                                        const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
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
                 "MultiParticleBrownianMotion.h\n";
#endif

  SetMolInBox(bPick);
  if (moleculeIndex.size() == 0) {
    std::cout << "Warning: MultiParticleBrownian move has no particles to move."
                 " Skipping...\n";
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
    return state;
  }

  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);
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

  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);

  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
  return state;
}

// To relax the system in NE_MTMC move
inline uint MultiParticleBrownian::PrepNEMTMC(const uint box, const uint midx,
                                              const uint kidx) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
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
    std::cout << "Warning: MultiParticleBrownian move has no particles to move."
                 " Skipping...\n";
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

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_MULTIPARTICLE_BM);
  return state;
}

inline uint MultiParticleBrownian::ChooseBox() {
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

inline uint MultiParticleBrownian::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_MULTIPARTICLE_BM);
  // Based on the reference force decide whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;

#ifdef GOMC_CUDA
  // This kernel will calculate translation/rotation amount + shifting/rotating
  if (moveType == mp::MPROTATE) {
    double r_max = moveSetRef.GetRMAX(bPick);
    BrownianMotionRotateParticlesGPU(
        cudaVars, moleculeIndex, molTorqueRef, newMolsPos, newCOMs, rt_k,
        boxDimRef.GetAxis(bPick), BETA, r_max, r123wrapper.GetStep(),
        r123wrapper.GetKeyValue(), r123wrapper.GetSeedValue(), bPick,
        isOrthogonal);
  } else {
    double t_max = moveSetRef.GetTMAX(bPick);
    BrownianMotionTranslateParticlesGPU(
        cudaVars, moleculeIndex, molForceRef, molForceRecRef, newMolsPos,
        newCOMs, rt_k, boxDimRef.GetAxis(bPick), BETA, t_max,
        r123wrapper.GetStep(), r123wrapper.GetKeyValue(),
        r123wrapper.GetSeedValue(), bPick, isOrthogonal);
  }
#else
  // Calculate trial translate and rotate
  // move particles according to force and torque and store them in the new pos
  CalculateTrialDistRot();
#endif

  // Do error checking and skip if there is an invalid transform amount.
  for (int m = 0; m < static_cast<int>(moleculeIndex.size()); ++m) {
    int molIndex = moleculeIndex[m];
    XYZ num = rt_k.Get(molIndex);
    if (moveType == mp::MPDISPLACE) { // displace
      // check for PBC error and bad initial configuration
      if (num > boxDimRef.GetHalfAxis(bPick)) {
        std::cout << "Warning: Trial Displacement exceeds half the box length "
                     "in the Multiparticle Brownian Motion move!\n";
        std::cout << "         Trial transformation vector: " << num << "\n";
        std::cout << "         Box Dimensions: " << boxDimRef.GetAxis(bPick);
        std::cout << "\n\n";
#ifndef NDEBUG
        std::cout << "Problem with molecule " << molIndex << std::endl;
#endif
        std::cout << "This might be due to a bad initial configuration, where "
                     "atoms of the molecules are too\n"
                  << "close to each other or overlap. Please equilibrate your "
                     "system using rigid body\n"
                  << "translation or rotation MC moves before using the "
                     "Multiparticle Brownian Motion move.\n\n";
        state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
        break;
      }
    }
    if (!std::isfinite(num.Length())) {
      std::cout << "Trial Displacement is not a finite number in Brownian"
                   " Motion Multiparticle move.\n";
      std::cout << "         Trial transform: " << num << "\n";
      state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
      break;
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_MULTIPARTICLE_BM);
  return state;
}

inline void MultiParticleBrownian::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_MULTIPARTICLE_BM);
  // Calculate the new force and energy and we will compare that to the
  // reference values in the Accept() function
  cellList.GridBox(boxDimRef, newMolsPos, molLookup, bPick);

  // back up cached Fourier term
  calcEwald->backupMolCache();

  // setup reciprocal vectors for new positions
  calcEwald->BoxReciprocalSums(bPick, newMolsPos);

  // calculate short range energy and force
  // this updates the Lennard-Jones and real electrostatic energy and
  // assigns the result to the new system potential
  sysPotNew = calcEnRef.BoxForce(sysPotRef, newMolsPos, atomForceRef,
                                 molForceNew, boxDimRef, moveType, bPick, true);

  // Calculate long range electrostatic energy for new positions
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

inline double MultiParticleBrownian::CalculateWRatio(XYZ const &lb_new,
                                                     XYZ const &lb_old,
                                                     XYZ const &k,
                                                     double max4) {
  double w_ratio = 0.0;
  XYZ old_var = lb_old - k;
  XYZ new_var = lb_new + k;

  // Note: we could factor max4 and multiply at the end, but
  //       for the move, where we translate and rotate all molecules,
  //       this method would not work. Hence, I did not factor it.
  //       It actually is w_ratio += -1.0* but we simplify it.
  w_ratio -= (new_var.LengthSq() / max4);
  // its actually is w_ratio -= -1.0* but we simplify it
  w_ratio += (old_var.LengthSq() / max4);

  return w_ratio;
}

inline double MultiParticleBrownian::GetCoeff() {
  // calculate (w_new->old/w_old->new) and return it.
  double w_ratio = 0.0;
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);
  double r_max4 = 4.0 * r_max;
  double t_max4 = 4.0 * t_max;

  if (moveType == mp::MPROTATE) { // rotate,
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(r_max, t_max, r_max4, t_max4) reduction(+:w_ratio)
#endif
    for (int m = 0; m < static_cast<int>(moleculeIndex.size()); ++m) {
      int molNumber = moleculeIndex[m];
      // rotate, bt_ = BETA * force * maxForce
      XYZ bt_old = molTorqueRef.Get(molNumber) * BETA * r_max;
      XYZ bt_new = molTorqueNew.Get(molNumber) * BETA * r_max;
      w_ratio += CalculateWRatio(bt_new, bt_old, rt_k.Get(molNumber), r_max4);
    }
  } else { // displace
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(r_max, t_max, r_max4, t_max4) reduction(+:w_ratio)
#endif
    for (int m = 0; m < static_cast<int>(moleculeIndex.size()); ++m) {
      int molNumber = moleculeIndex[m];
      // bf_ = BETA * torque * maxTorque
      XYZ bf_old =
          (molForceRef.Get(molNumber) + molForceRecRef.Get(molNumber)) * BETA *
          t_max;
      XYZ bf_new =
          (molForceNew.Get(molNumber) + molForceRecNew.Get(molNumber)) * BETA *
          t_max;
      w_ratio += CalculateWRatio(bf_new, bf_old, rt_k.Get(molNumber), t_max4);
    }
  }

  return w_ratio;
}

inline void MultiParticleBrownian::Accept(const uint rejectState,
                                          const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_MULTIPARTICLE_BM);
  // Compare the energies of reference and trial positions to decide whether to
  // accept or reject the move
  double MPCoeff = GetCoeff();
  double accept =
      exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()) + MPCoeff);
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if (result) {
    sysPotRef = sysPotNew;
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    // Update reciprocal value
    calcEwald->UpdateRecip(bPick);
    // Update the velocity in box
    velocity.UpdateBoxVelocity(bPick);
  } else {
    cellList.GridBox(boxDimRef, coordCurrRef, molLookup, bPick);
    calcEwald->exgMolCache();
  }

  moveSetRef.UpdateMoveSettingMultiParticle(bPick, result, moveType);
  moveSetRef.Update(mv::MULTIPARTICLE_BM, result, bPick);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_MULTIPARTICLE_BM);
}

inline XYZ MultiParticleBrownian::CalcRandomTransform(XYZ const &lb,
                                                      double const max,
                                                      int molIndex) {
  XYZ lbmax = lb * max;
  XYZ num, randnums;
  // variance is 2A according to the paper, so stdDev is sqrt(variance)
  double stdDev = sqrt(2.0 * max);

  randnums = r123wrapper.GetGaussianCoords(molIndex, 0.0, stdDev);
  num.x = lbmax.x + randnums.x;
  num.y = lbmax.y + randnums.y;
  num.z = lbmax.z + randnums.z;

  return num;
}

inline void MultiParticleBrownian::CalculateTrialDistRot() {
  double r_max = moveSetRef.GetRMAX(bPick);
  double t_max = moveSetRef.GetTMAX(bPick);
  double *x = rt_k.x;
  double *y = rt_k.y;
  double *z = rt_k.z;
  if (moveType == mp::MPROTATE) { // rotate
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(r_max, x, y, z)
#endif
    for (int m = 0; m < static_cast<int>(moleculeIndex.size()); ++m) {
      int molIndex = moleculeIndex[m];
      XYZ bt = molTorqueRef.Get(molIndex) * BETA;
      XYZ val = CalcRandomTransform(bt, r_max, molIndex);
      x[molIndex] = val.x;
      y[molIndex] = val.y;
      z[molIndex] = val.z;
      RotateForceBiased(molIndex);
    }
  } else { // displace
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(t_max, x, y, z)
#endif
    for (int m = 0; m < static_cast<int>(moleculeIndex.size()); ++m) {
      int molIndex = moleculeIndex[m];
      XYZ bf =
          (molForceRef.Get(molIndex) + molForceRecRef.Get(molIndex)) * BETA;
      XYZ val = CalcRandomTransform(bf, t_max, molIndex);
      x[molIndex] = val.x;
      y[molIndex] = val.y;
      z[molIndex] = val.z;
      TranslateForceBiased(molIndex);
    }
  }
}

inline void MultiParticleBrownian::RotateForceBiased(int molIndex) {
  XYZ rot = rt_k.Get(molIndex);
  double rotLen = rot.Length();
  RotationMatrix matrix;

  matrix = RotationMatrix::FromAxisAngle(rotLen, rot.Normalize());

  XYZ center = comCurrRef.Get(molIndex);
  XYZ CoMoffset = center - matrix.Apply(center);
  uint start, stop, len;
  molRef.GetRange(start, stop, len, molIndex);

  // Copy the range into temporary array
  XYZArray temp(len);
  newMolsPos.CopyRange(temp, start, 0, len);
  boxDimRef.UnwrapPBC(temp, bPick, center);

  // Computing the rotation this way has less round off than other methods
  // Rotate each atom in the molecule
  for (int p = 0; p < static_cast<int>(len); ++p) {
    temp.Set(p, matrix.Apply(temp[p]));
    temp.Add(p, CoMoffset);
  }
  boxDimRef.WrapPBC(temp, bPick);
  // Copy back the result
  temp.CopyRange(newMolsPos, 0, start, len);
}

inline void MultiParticleBrownian::TranslateForceBiased(int molIndex) {
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

#endif
