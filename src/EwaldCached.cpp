/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "EwaldCached.h"

#include "GOMCEventsProfile.h"
#include "StaticVals.h"

using namespace geom;

template <class T> void SafeDeleteArray(T *&pVal) {
  delete[] pVal;
  pVal = NULL;
}

EwaldCached::EwaldCached(StaticVals &stat, System &sys)
    : Ewald(stat, sys)
#if ENSEMBLE == GEMC
      ,
      GEMC_KIND(stat.kindOfGEMC)
#endif
{
}

EwaldCached::~EwaldCached() {
#ifdef _OPENMP
#pragma omp parallel for default(none)
#endif
  for (int i = 0; i < (int)mols.count; i++) {
    SafeDeleteArray(cosMolRef[i]);
    SafeDeleteArray(sinMolRef[i]);
    SafeDeleteArray(cosMolBoxRecip[i]);
    SafeDeleteArray(sinMolBoxRecip[i]);
  }

  SafeDeleteArray(cosMolRef);
  SafeDeleteArray(sinMolRef);
  SafeDeleteArray(cosMolBoxRecip);
  SafeDeleteArray(sinMolBoxRecip);

  SafeDeleteArray(cosMolRestore);
  SafeDeleteArray(sinMolRestore);
}

void EwaldCached::Init() {
  for (uint m = 0; m < mols.count; ++m) {
    const MoleculeKind &molKind = mols.GetKind(m);
    for (uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleKind.push_back(molKind.AtomKind(a));
      particleMol.push_back(m);
      particleCharge.push_back(molKind.AtomCharge(a));
      if (std::abs(molKind.AtomCharge(a)) < 1e-9) {
        particleHasNoCharge.push_back(true);
      } else {
        particleHasNoCharge.push_back(false);
      }
    }
  }

  AllocMem();
  // initialize K vectors and reciprocal terms
  UpdateVectorsAndRecipTerms(true);
}

void EwaldCached::AllocMem() {
  // Use the base class member function to allocate base class data members
  Ewald::AllocMem();

  cosMolRestore = new double[imageTotal];
  sinMolRestore = new double[imageTotal];
  cosMolRef = new double *[mols.count];
  sinMolRef = new double *[mols.count];
  cosMolBoxRecip = new double *[mols.count];
  sinMolBoxRecip = new double *[mols.count];

#ifdef _OPENMP
#pragma omp parallel for default(none)
#endif
  for (int i = 0; i < (int)mols.count; i++) {
    cosMolRef[i] = new double[imageTotal];
    sinMolRef[i] = new double[imageTotal];
    cosMolBoxRecip[i] = new double[imageTotal];
    sinMolBoxRecip[i] = new double[imageTotal];
  }
}

// calculate reciprocal terms for a box. Should be called only at
// the start of the simulation to initialize the settings and when
// testing a change in box dimensions, such as a volume transfer.
void EwaldCached::BoxReciprocalSetup(uint box, XYZArray const &molCoords) {
  if (box < BOXES_WITH_U_NB) {
    GOMC_EVENT_START(1, GomcProfileEvent::RECIP_BOX_SETUP);
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

#ifdef _OPENMP
#pragma omp parallel sections default(none) shared(box)
    {
#pragma omp section
      std::memset(sumRnew[box], 0.0, sizeof(double) * imageSize[box]);
#pragma omp section
      std::memset(sumInew[box], 0.0, sizeof(double) * imageSize[box]);
    }
#else
    std::memset(sumRnew[box], 0.0, sizeof(double) * imageSize[box]);
    std::memset(sumInew[box], 0.0, sizeof(double) * imageSize[box]);
#endif

    while (thisMol != end) {
      MoleculeKind const &thisKind = mols.GetKind(*thisMol);
      double lambdaCoef = GetLambdaCoef(*thisMol, box);
      uint startAtom = mols.MolStart(*thisMol);

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(box, lambdaCoef, molCoords, startAtom, thisKind, thisMol)
#endif
      for (int i = 0; i < (int)imageSize[box]; i++) {
        cosMolRef[*thisMol][i] = 0.0;
        sinMolRef[*thisMol][i] = 0.0;

        for (uint j = 0; j < thisKind.NumAtoms(); j++) {
          if (particleHasNoCharge[startAtom + j]) {
            continue;
          }
          double dotProduct = Dot(mols.MolStart(*thisMol) + j, kx[box][i],
                                  ky[box][i], kz[box][i], molCoords);
          // cache the Cos and sin term with lambda = 1
          cosMolRef[*thisMol][i] += (thisKind.AtomCharge(j) * cos(dotProduct));
          sinMolRef[*thisMol][i] += (thisKind.AtomCharge(j) * sin(dotProduct));
        }
        // store the summation with system lambda
        sumRnew[box][i] += (lambdaCoef * cosMolRef[*thisMol][i]);
        sumInew[box][i] += (lambdaCoef * sinMolRef[*thisMol][i]);
      }

      thisMol++;
    }
    GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_BOX_SETUP);
  }
}

// Calculate reciprocal terms for a box, when an updated value is needed
// because the number and location of molecules could have changed since
// the volume was set. Examples include MultiParticle and Molecule Exchange
// moves. For these calls, we need to use the Reference settings, since
// these hold the information for the current box dimensions.
void EwaldCached::BoxReciprocalSums(uint box, XYZArray const &molCoords) {
  if (box < BOXES_WITH_U_NB) {
    GOMC_EVENT_START(1, GomcProfileEvent::RECIP_BOX_SETUP);
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

#ifdef _OPENMP
#pragma omp parallel sections default(none) shared(box)
    {
#pragma omp section
      std::memset(sumRnew[box], 0.0, sizeof(double) * imageSizeRef[box]);
#pragma omp section
      std::memset(sumInew[box], 0.0, sizeof(double) * imageSizeRef[box]);
    }
#else
    std::memset(sumRnew[box], 0.0, sizeof(double) * imageSizeRef[box]);
    std::memset(sumInew[box], 0.0, sizeof(double) * imageSizeRef[box]);
#endif

    while (thisMol != end) {
      MoleculeKind const &thisKind = mols.GetKind(*thisMol);
      double lambdaCoef = GetLambdaCoef(*thisMol, box);
      uint startAtom = mols.MolStart(*thisMol);

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(box, lambdaCoef, molCoords, startAtom, thisKind, thisMol)
#endif
      for (int i = 0; i < (int)imageSizeRef[box]; i++) {
        cosMolRef[*thisMol][i] = 0.0;
        sinMolRef[*thisMol][i] = 0.0;

        for (uint j = 0; j < thisKind.NumAtoms(); j++) {
          if (particleHasNoCharge[startAtom + j]) {
            continue;
          }
          double dotProduct = Dot(mols.MolStart(*thisMol) + j, kxRef[box][i],
                                  kyRef[box][i], kzRef[box][i], molCoords);
          // cache the Cos and sin term with lambda = 1
          cosMolRef[*thisMol][i] += (thisKind.AtomCharge(j) * cos(dotProduct));
          sinMolRef[*thisMol][i] += (thisKind.AtomCharge(j) * sin(dotProduct));
        }
        // store the summation with system lambda
        sumRnew[box][i] += (lambdaCoef * cosMolRef[*thisMol][i]);
        sumInew[box][i] += (lambdaCoef * sinMolRef[*thisMol][i]);
      }

      thisMol++;
    }
    GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_BOX_SETUP);
  }
}

// Calculate reciprocal energy for a box. This function is called for two
// reasons:
// 1. During a volume move, we call this function to total the sums for the new
//    volume. For these calls, need to total the sums with the new volume
//    settings.
// 2. During a Molecule Exchange or MultiParticle move, we need to recompute
// these
//    sums for the current volume, since the number and location of molecules
//    could have changed since the volume was set. For these calls, we need to
//    use the Reference settings, since these hold the information for the
//    current box dimensions. Also called at the start of the simulation, after
//    the Reference volume parameters have been set.
double EwaldCached::BoxReciprocal(uint box, bool isNewVolume) const {
  double energyRecip = 0.0;

  if (box < BOXES_WITH_U_NB) {
    // GOMC_EVENT_START(1, GomcProfileEvent::RECIP_BOX_ENERGY);
    double *prefactPtr;
    int imageSzVal;
    if (isNewVolume) {
      prefactPtr = prefact[box];
      imageSzVal = static_cast<int>(imageSize[box]);
    } else {
      prefactPtr = prefactRef[box];
      imageSzVal = static_cast<int>(imageSizeRef[box]);
    }
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(box, imageSzVal, prefactPtr)     \
    reduction(+ : energyRecip)
#endif
    for (int i = 0; i < imageSzVal; i++) {
      energyRecip += ((sumRnew[box][i] * sumRnew[box][i] +
                       sumInew[box][i] * sumInew[box][i]) *
                      prefactPtr[i]);
    }
    // GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_BOX_ENERGY);
  }

  return energyRecip;
}

// calculate reciprocal term for displacement and rotation move
double EwaldCached::MolReciprocal(XYZArray const &molCoords,
                                  const uint molIndex, const uint box) {
  double energyRecipNew = 0.0;

  if (box < BOXES_WITH_U_NB) {
    GOMC_EVENT_START(1, GomcProfileEvent::RECIP_MOL_ENERGY);
    MoleculeKind const &thisKind = mols.GetKind(molIndex);
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);
    double lambdaCoef = GetLambdaCoef(molIndex, box);

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(lambdaCoef, length, molCoords, startAtom, thisKind)                 \
    reduction(+ : energyRecipNew) firstprivate(box, molIndex)
#endif
    for (int i = 0; i < (int)imageSizeRef[box]; i++) {
      double sumRealNew = 0.0;
      double sumImaginaryNew = 0.0;
      double sumRealOld = cosMolRef[molIndex][i];
      double sumImaginaryOld = sinMolRef[molIndex][i];
      cosMolRestore[i] = cosMolRef[molIndex][i];
      sinMolRestore[i] = sinMolRef[molIndex][i];

      for (uint p = 0; p < length; ++p) {
        if (particleHasNoCharge[startAtom + p]) {
          continue;
        }
        double dotProductNew =
            Dot(p, kxRef[box][i], kyRef[box][i], kzRef[box][i], molCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));
      }

      sumRnew[box][i] =
          sumRref[box][i] + lambdaCoef * (sumRealNew - sumRealOld);
      sumInew[box][i] =
          sumIref[box][i] + lambdaCoef * (sumImaginaryNew - sumImaginaryOld);
      cosMolRef[molIndex][i] = sumRealNew;
      sinMolRef[molIndex][i] = sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] +
                         sumInew[box][i] * sumInew[box][i]) *
                        prefactRef[box][i];
    }
    GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_MOL_ENERGY);
  }

  return energyRecipNew - sysPotRef.boxEnergy[box].recip;
}

// calculate reciprocal term in destination box for swap move
// No need to scale the charge with lambda, since this function will not be
// called in free energy of NeMTMC
double EwaldCached::SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                                  const int molIndex) {
  GOMC_EVENT_START(1, GomcProfileEvent::RECIP_SWAP_ENERGY);
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

#ifdef _OPENMP
#pragma omp parallel sections default(none) firstprivate(molIndex)
  {
#pragma omp section
    std::memcpy(cosMolRestore, cosMolRef[molIndex],
                sizeof(double) * imageTotal);
#pragma omp section
    std::memcpy(sinMolRestore, sinMolRef[molIndex],
                sizeof(double) * imageTotal);
  }
#else
  std::memcpy(cosMolRestore, cosMolRef[molIndex], sizeof(double) * imageTotal);
  std::memcpy(sinMolRestore, sinMolRef[molIndex], sizeof(double) * imageTotal);
#endif

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const &thisKind = newMol.GetKind();
    XYZArray molCoords = newMol.GetCoords();
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(length, molCoords, startAtom, thisKind) firstprivate(box, molIndex) \
    reduction(+ : energyRecipNew)
#endif
    for (int i = 0; i < (int)imageSizeRef[box]; i++) {
      cosMolRef[molIndex][i] = 0.0;
      sinMolRef[molIndex][i] = 0.0;

      for (uint p = 0; p < length; ++p) {
        if (particleHasNoCharge[startAtom + p]) {
          continue;
        }
        double dotProductNew =
            Dot(p, kxRef[box][i], kyRef[box][i], kzRef[box][i], molCoords);
        cosMolRef[molIndex][i] += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sinMolRef[molIndex][i] += (thisKind.AtomCharge(p) * sin(dotProductNew));
      }

      // sumRealNew;
      sumRnew[box][i] = sumRref[box][i] + cosMolRef[molIndex][i];
      // sumImaginaryNew;
      sumInew[box][i] = sumIref[box][i] + sinMolRef[molIndex][i];

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] +
                         sumInew[box][i] * sumInew[box][i]) *
                        prefactRef[box][i];
    }

    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_SWAP_ENERGY);
  return energyRecipNew - energyRecipOld;
}

// calculate reciprocal term in source box for swap move
// No need to scale the charge with lambda, since this function will not be
// called in free energy of NeMTMC
double EwaldCached::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                    const uint box, const int molIndex) {
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

  if (box < BOXES_WITH_U_NB) {
    GOMC_EVENT_START(1, GomcProfileEvent::RECIP_SWAP_ENERGY);
#ifdef _OPENMP
#pragma omp parallel for default(none) firstprivate(box)                       \
    reduction(+ : energyRecipNew)
#endif
    for (int i = 0; i < (int)imageSizeRef[box]; i++) {
      sumRnew[box][i] = sumRref[box][i] - cosMolRestore[i];
      sumInew[box][i] = sumIref[box][i] - sinMolRestore[i];

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] +
                         sumInew[box][i] * sumInew[box][i]) *
                        prefactRef[box][i];
    }

    energyRecipOld = sysPotRef.boxEnergy[box].recip;
    GOMC_EVENT_STOP(1, GomcProfileEvent::RECIP_SWAP_ENERGY);
  }
  return energyRecipNew - energyRecipOld;
}

// calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double
EwaldCached::MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                                   const std::vector<cbmc::TrialMol> &oldMol,
                                   const std::vector<uint> &molIndexNew,
                                   const std::vector<uint> &molIndexOld) {
  // This function should not be called in IDExchange move
  std::cout << "Error: Cached Fourier method cannot be used while "
            << "performing Molecule Exchange move!" << std::endl;
  exit(EXIT_FAILURE);
  return 0.0;
}

// calculate reciprocal term for lambdaNew and Old with same coordinates
double EwaldCached::ChangeLambdaRecip(XYZArray const &molCoords,
                                      const double lambdaOld,
                                      const double lambdaNew,
                                      const uint molIndex, const uint box) {
  // This function should not be called in NeMTMC move
  std::cout << "Error: Cached Fourier method cannot be used while "
            << "performing non-equilibrium Mol-Transfer MC move (NeMTMC)!"
            << std::endl;
  exit(EXIT_FAILURE);
  return 0.0;
}

// calculate reciprocal term for lambdaNew and Old with same coordinates
void EwaldCached::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                              const std::vector<double> &lambda_Coul,
                              const uint iState, const uint molIndex,
                              const uint box) const {
  // Need to implement GPU
  uint lambdaSize = lambda_Coul.size();
  double *energyRecip = new double[lambdaSize];
  std::fill_n(energyRecip, lambdaSize, 0.0);

#if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
#pragma omp parallel for default(none) shared(lambda_Coul, lambdaSize)         \
    reduction(+ : energyRecip[ : lambdaSize])                                  \
    firstprivate(box, iState, molIndex)
#endif
  for (uint i = 0; i < imageSizeRef[box]; i++) {
    for (uint s = 0; s < lambdaSize; s++) {
      // Calculate the energy of other state
      double coefDiff = sqrt(lambda_Coul[s]) - sqrt(lambda_Coul[iState]);
      energyRecip[s] +=
          prefactRef[box][i] *
          ((sumRref[box][i] + coefDiff * cosMolRef[molIndex][i]) *
               (sumRref[box][i] + coefDiff * cosMolRef[molIndex][i]) +
           (sumIref[box][i] + coefDiff * sinMolRef[molIndex][i]) *
               (sumIref[box][i] + coefDiff * sinMolRef[molIndex][i]));
    }
  }

  double energyRecipOld = sysPotRef.boxEnergy[box].recip;
  for (uint s = 0; s < lambdaSize; s++) {
    energyDiff[s].recip = energyRecip[s] - energyRecipOld;
  }
  // Calculate du/dl of Reciprocal for current state
  // energy difference E(lambda =1) - E(lambda = 0)
  dUdL_Coul.recip += energyDiff[lambdaSize - 1].recip - energyDiff[0].recip;
  delete[] energyRecip;
}

// restore cosMol and sinMol
void EwaldCached::RestoreMol(int molIndex) {
  double *tempCos, *tempSin;
  tempCos = cosMolRef[molIndex];
  tempSin = sinMolRef[molIndex];
  cosMolRef[molIndex] = cosMolRestore;
  sinMolRef[molIndex] = sinMolRestore;
  cosMolRestore = tempCos;
  sinMolRestore = tempSin;
}

// restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void EwaldCached::exgMolCache() {
  double **tempCos, **tempSin;
  tempCos = cosMolRef;
  tempSin = sinMolRef;
  cosMolRef = cosMolBoxRecip;
  sinMolRef = sinMolBoxRecip;
  cosMolBoxRecip = tempCos;
  sinMolBoxRecip = tempSin;
}

// backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void EwaldCached::backupMolCache() {
#if ENSEMBLE == NPT || ENSEMBLE == NVT
  exgMolCache();
#else
#ifdef _OPENMP
#pragma omp parallel for default(none)
#endif
  for (int m = 0; m < (int)mols.count; m++) {
    std::memcpy(cosMolBoxRecip[m], cosMolRef[m], sizeof(double) * imageTotal);
    std::memcpy(sinMolBoxRecip[m], sinMolRef[m], sizeof(double) * imageTotal);
  }
#endif
}
