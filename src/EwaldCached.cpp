/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "EwaldCached.h"
#include "StaticVals.h"

using namespace geom;


template< class T > void SafeDeleteArray( T*& pVal )
{
  delete[] pVal;
  pVal = NULL;
}

EwaldCached::EwaldCached(StaticVals & stat, System & sys) : Ewald(stat, sys)
#if ENSEMBLE == GEMC
  , GEMC_KIND(stat.kindOfGEMC)
#endif
{}


EwaldCached::~EwaldCached()
{
#ifdef _OPENMP
  #pragma omp parallel for default(none)
#endif
  for (int i = 0; i < mols.count; i++) {
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

void EwaldCached::Init()
{
  for(uint m = 0; m < mols.count; ++m) {
    const MoleculeKind& molKind = mols.GetKind(m);
    for(uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleKind.push_back(molKind.AtomKind(a));
      particleMol.push_back(m);
      particleCharge.push_back(molKind.AtomCharge(a));
      if(std::abs(molKind.AtomCharge(a)) < 1e-9) {
        particleHasNoCharge.push_back(true);
      } else {
        particleHasNoCharge.push_back(false);
      }
    }
  }

  AllocMem();
  //initialize K vectors and reciprocate terms
  UpdateVectorsAndRecipTerms(true);
}

void EwaldCached::AllocMem()
{
  //Use the base class member function to allocate base class data members
  Ewald::AllocMem();

  cosMolRestore = new double[imageTotal];
  sinMolRestore = new double[imageTotal];
  cosMolRef = new double*[mols.count];
  sinMolRef = new double*[mols.count];
  cosMolBoxRecip = new double*[mols.count];
  sinMolBoxRecip = new double*[mols.count];

#ifdef _OPENMP
  #pragma omp parallel for default(none)
#endif
  for (int i = 0; i < mols.count; i++) {
    cosMolRef[i] = new double[imageTotal];
    sinMolRef[i] = new double[imageTotal];
    cosMolBoxRecip[i] = new double[imageTotal];
    sinMolBoxRecip[i] = new double[imageTotal];
  }
}

//calculate reciprocal term for a box
void EwaldCached::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{

  if (box < BOXES_WITH_U_NB) {
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);

#ifdef _OPENMP
    #pragma omp parallel default(none) shared(box)
#endif
    {
      std::memset(sumRnew[box], 0.0, sizeof(double) * imageSize[box]);
      std::memset(sumInew[box], 0.0, sizeof(double) * imageSize[box]);
    }

    while (thisMol != end) {
      MoleculeKind const& thisKind = mols.GetKind(*thisMol);
      double lambdaCoef = GetLambdaCoef(*thisMol, box);
      uint startAtom = mols.MolStart(*thisMol);

#ifdef _OPENMP
      #pragma omp parallel for default(none) shared(box, lambdaCoef, molCoords, startAtom, thisKind, thisMol)
#endif
      for (int i = 0; i < imageSize[box]; i++) {
        cosMolRef[*thisMol][i] = 0.0;
        sinMolRef[*thisMol][i] = 0.0;

        for (uint j = 0; j < thisKind.NumAtoms(); j++) {
          if(particleHasNoCharge[startAtom + j]) {
            continue;
          }
          double dotProduct = Dot(mols.MolStart(*thisMol) + j,
                                  kx[box][i], ky[box][i],
                                  kz[box][i], molCoords);
          //cache the Cos and sin term with lambda = 1
          cosMolRef[*thisMol][i] += (thisKind.AtomCharge(j) *
                                     cos(dotProduct));
          sinMolRef[*thisMol][i] += (thisKind.AtomCharge(j) *
                                     sin(dotProduct));
        }
        //store the summation with system lambda
        sumRnew[box][i] += (lambdaCoef * cosMolRef[*thisMol][i]);
        sumInew[box][i] += (lambdaCoef * sinMolRef[*thisMol][i]);
      }

      thisMol++;
    }
  }
}


//calculate reciprocal term for a box
double EwaldCached::BoxReciprocal(uint box) const
{
  double energyRecip = 0.0;

  if (box < BOXES_WITH_U_NB) {
#ifdef _OPENMP
    #pragma omp parallel for default(none) shared(box) reduction(+:energyRecip)
#endif
    for (int i = 0; i < imageSize[box]; i++) {
      energyRecip += (( sumRnew[box][i] * sumRnew[box][i] +
                        sumInew[box][i] * sumInew[box][i]) *
                      prefact[box][i]);
    }
  }

  return energyRecip;
}

//calculate reciprocal term for displacement and rotation move
double EwaldCached::MolReciprocal(XYZArray const& molCoords,
                                  const uint molIndex,
                                  const uint box)
{
  double energyRecipNew = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const& thisKind = mols.GetKind(molIndex);
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);
    double lambdaCoef = GetLambdaCoef(molIndex, box);

#ifdef _OPENMP
#if GCC_VERSION >= 90000
    #pragma omp parallel for default(none) shared(lambdaCoef, length, molCoords, \
    startAtom, thisKind, box, molIndex) \
reduction(+:energyRecipNew)
#else
    #pragma omp parallel for default(none) shared(lambdaCoef, length, molCoords, \
    startAtom, thisKind) \
reduction(+:energyRecipNew)
#endif
#endif
    for (int i = 0; i < imageSizeRef[box]; i++) {
      double sumRealNew = 0.0;
      double sumImaginaryNew = 0.0;
      double sumRealOld = cosMolRef[molIndex][i];
      double sumImaginaryOld = sinMolRef[molIndex][i];
      cosMolRestore[i] = cosMolRef[molIndex][i];
      sinMolRestore[i] = sinMolRef[molIndex][i];

      for (uint p = 0; p < length; ++p) {
        if(particleHasNoCharge[startAtom + p]) {
          continue;
        }
        double dotProductNew = Dot(p, kxRef[box][i],
                                   kyRef[box][i], kzRef[box][i],
                                   molCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));
      }

      sumRnew[box][i] = sumRref[box][i] + lambdaCoef *
                        (sumRealNew - sumRealOld);
      sumInew[box][i] = sumIref[box][i] + lambdaCoef *
                        (sumImaginaryNew - sumImaginaryOld);
      cosMolRef[molIndex][i] = sumRealNew;
      sinMolRef[molIndex][i] = sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
  }

  return energyRecipNew - sysPotRef.boxEnergy[box].recip;
}

//calculate reciprocal term in destination box for swap move
//No need to scale the charge with lambda, since this function will not be
// called in free energy of CFCMC
double EwaldCached::SwapDestRecip(const cbmc::TrialMol &newMol,
                                  const uint box,
                                  const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

#ifdef _OPENMP
#if GCC_VERSION >= 90000
  #pragma omp parallel default(none) shared(molIndex)
#else
  #pragma omp parallel default(none)
#endif
#endif
  {
    std::memcpy(cosMolRestore, cosMolRef[molIndex], sizeof(double)*imageTotal);
    std::memcpy(sinMolRestore, sinMolRef[molIndex], sizeof(double)*imageTotal);
  }

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const& thisKind = newMol.GetKind();
    XYZArray molCoords = newMol.GetCoords();
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);

#ifdef _OPENMP
#if GCC_VERSION >= 90000
    #pragma omp parallel for default(none) shared(length, molCoords, startAtom, \
    thisKind, box, molIndex) \
reduction(+:energyRecipNew)
#else
    #pragma omp parallel for default(none) shared(length, molCoords, startAtom, \
    thisKind) \
reduction(+:energyRecipNew)
#endif
#endif
    for (int i = 0; i < imageSizeRef[box]; i++) {
      cosMolRef[molIndex][i] = 0.0;
      sinMolRef[molIndex][i] = 0.0;

      for (uint p = 0; p < length; ++p) {
        if(particleHasNoCharge[startAtom + p]) {
          continue;
        }
        double dotProductNew = Dot(p, kxRef[box][i],
                                   kyRef[box][i], kzRef[box][i],
                                   molCoords);
        cosMolRef[molIndex][i] += (thisKind.AtomCharge(p) *
                                   cos(dotProductNew));
        sinMolRef[molIndex][i] += (thisKind.AtomCharge(p) *
                                   sin(dotProductNew));
      }

      //sumRealNew;
      sumRnew[box][i] = sumRref[box][i] + cosMolRef[molIndex][i];
      //sumImaginaryNew;
      sumInew[box][i] = sumIref[box][i] + sinMolRef[molIndex][i];

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }

    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }

  return energyRecipNew - energyRecipOld;
}


//calculate reciprocal term in source box for swap move
//No need to scale the charge with lambda, since this function will not be
// called in free energy of CFCMC
double EwaldCached::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                    const uint box, const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

  if (box < BOXES_WITH_U_NB) {
#ifdef _OPENMP
#if GCC_VERSION >= 90000
    #pragma omp parallel for default(none) shared(box) reduction(+:energyRecipNew)
#else
    #pragma omp parallel for default(none) reduction(+:energyRecipNew)
#endif
#endif
    for (int i = 0; i < imageSizeRef[box]; i++) {
      sumRnew[box][i] = sumRref[box][i] - cosMolRestore[i];
      sumInew[box][i] = sumIref[box][i] - sinMolRestore[i];

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }

    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }
  return energyRecipNew - energyRecipOld;
}

//calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double EwaldCached::SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
                              const std::vector<cbmc::TrialMol> &oldMol,
                              const std::vector<uint> molIndexNew,
                              const std::vector<uint> molIndexOld)
{
  //This function should not be called in IDExchange move
  std::cout << "Error: Cached Fourier method cannot be used while " <<
            "performing Molecule Exchange move!" << std::endl;
  exit(EXIT_FAILURE);
  return 0.0;
}

//calculate reciprocal term for lambdaNew and Old with same coordinates
double EwaldCached::CFCMCRecip(XYZArray const& molCoords, const double lambdaOld,
                               const double lambdaNew, const uint molIndex,
                               const uint box)
{
  //This function should not be called in CFCMC move
  std::cout << "Error: Cached Fourier method cannot be used while " <<
            "performing CFCMC move!" << std::endl;
  exit(EXIT_FAILURE);
  return 0.0;
}

//calculate reciprocal term for lambdaNew and Old with same coordinates
void EwaldCached::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                              const std::vector<double> &lambda_Coul,
                              const uint iState, const uint molIndex,
                              const uint box) const
{
  //Need to implement GPU
  uint lambdaSize = lambda_Coul.size();
  double *energyRecip = new double [lambdaSize];
  std::fill_n(energyRecip, lambdaSize, 0.0);

#if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
#if GCC_VERSION >= 90000
  #pragma omp parallel for default(none) shared(lambda_Coul, lambdaSize, box, \
  iState, molIndex) \
reduction(+:energyRecip[:lambdaSize])
#else
  #pragma omp parallel for default(none) shared(lambda_Coul, lambdaSize) \
  reduction(+:energyRecip[:lambdaSize])
#endif
#endif
  for (uint i = 0; i < imageSizeRef[box]; i++) {
    for(uint s = 0; s < lambdaSize; s++) {
      //Calculate the energy of other state
      double coefDiff = sqrt(lambda_Coul[s]) - sqrt(lambda_Coul[iState]);
      energyRecip[s] += prefactRef[box][i] *
                        ((sumRref[box][i] + coefDiff * cosMolRef[molIndex][i]) *
                         (sumRref[box][i] + coefDiff * cosMolRef[molIndex][i]) +
                         (sumIref[box][i] + coefDiff * sinMolRef[molIndex][i]) *
                         (sumIref[box][i] + coefDiff * sinMolRef[molIndex][i]));
    }
  }

  double energyRecipOld = sysPotRef.boxEnergy[box].recip;
  for(uint s = 0; s < lambdaSize; s++) {
  energyDiff[s].recip = energyRecip[s] - energyRecipOld;
  }
  //Calculate du/dl of Reciprocal for current state
  //energy difference E(lambda =1) - E(lambda = 0)
  dUdL_Coul.recip += energyDiff[lambdaSize - 1].recip - energyDiff[0].recip;
  delete [] energyRecip;
}

//restore cosMol and sinMol
void EwaldCached::RestoreMol(int molIndex)
{
  double *tempCos, *tempSin;
  tempCos = cosMolRef[molIndex];
  tempSin = sinMolRef[molIndex];
  cosMolRef[molIndex] = cosMolRestore;
  sinMolRef[molIndex] = sinMolRestore;
  cosMolRestore = tempCos;
  sinMolRestore = tempSin;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void EwaldCached::exgMolCache()
{
  double **tempCos, **tempSin;
  tempCos = cosMolRef;
  tempSin = sinMolRef;
  cosMolRef = cosMolBoxRecip;
  sinMolRef = sinMolBoxRecip;
  cosMolBoxRecip = tempCos;
  sinMolBoxRecip = tempSin;
}

//backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void EwaldCached::backupMolCache()
{
#if ENSEMBLE == NPT || ENSEMBLE == NVT
  exgMolCache();
#else
#ifdef _OPENMP
  #pragma omp parallel for default(none)
#endif
  for(int m = 0; m < mols.count; m++) {
    std::memcpy(cosMolBoxRecip[m], cosMolRef[m], sizeof(double) * imageTotal);
    std::memcpy(sinMolBoxRecip[m], sinMolRef[m], sizeof(double) * imageTotal);
  }
#endif
}
