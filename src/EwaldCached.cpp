/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "EwaldCached.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"             //
#include "MoleculeLookup.h"
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#endif

//
//
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald
//
//    Developed by Y. Li and Mohammad S. Barhaghi
//
//

using namespace geom;

EwaldCached::EwaldCached(StaticVals & stat, System & sys) : Ewald(stat, sys) { }

EwaldCached::~EwaldCached()
{
  if(ewald) {
    for(int i = 0; i < mols.count; i++) {
      //when cached option is choosen
      if (cosMolRef[i] != NULL) {
        delete[] cosMolRef[i];
        delete[] sinMolRef[i];
        delete[] cosMolBoxRecip[i];
        delete[] sinMolBoxRecip[i];
      }
    }

    if (kx != NULL) {
      //when cached option is choosen
      if (cosMolRestore != NULL) {
        delete[] cosMolRestore;
        delete[] sinMolRestore;
      }
      //when cached option is choosen
      if (cosMolRef != NULL) {
        delete[] cosMolRef;
        delete[] sinMolRef;
        delete[] cosMolBoxRecip;
        delete[] sinMolBoxRecip;
      }
    }
  }
}

void EwaldCached::Init()
{
  for(uint m = 0; m < mols.count; ++m) {
    const MoleculeKind& molKind = mols.GetKind(m);
    for(uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleKind.push_back(molKind.AtomKind(a));
      particleMol.push_back(m);
      particleCharge.push_back(molKind.AtomCharge(a));
    }
  }

  electrostatic = forcefield.electrostatic;
  ewald = forcefield.ewald;
  alpha = forcefield.alpha;
  recip_rcut = forcefield.recip_rcut;
  recip_rcut_Sq = recip_rcut * recip_rcut;
  AllocMem();
  //initialize K vectors and reciprocate terms
  UpdateVectorsAndRecipTerms();
}

void EwaldCached::AllocMem()
{
  cosMolRef = new double*[mols.count];
  sinMolRef = new double*[mols.count];
  cosMolBoxRecip = new double*[mols.count];
  sinMolBoxRecip = new double*[mols.count];

  //25% larger than original box size, reserved for image size change
  imageTotal = findLargeImage();

  cosMolRestore = new double[imageTotal];
  sinMolRestore = new double[imageTotal];

  int i;
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < mols.count; i++) {
    cosMolRef[i] = new double[imageTotal];
    sinMolRef[i] = new double[imageTotal];
    cosMolBoxRecip[i] = new double[imageTotal];
    sinMolBoxRecip[i] = new double[imageTotal];
  }

}

//calculate reciprocate term for a box
void EwaldCached::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
  uint j, m;
  int i;
  double dotProduct = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);

#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
      std::memset(sumRnew[box], 0.0, sizeof(double) * imageSize[box]);
      std::memset(sumInew[box], 0.0, sizeof(double) * imageSize[box]);
    }

    while (thisMol != end) {
      MoleculeKind const& thisKind = mols.GetKind(*thisMol);

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(i, j, dotProduct)
#endif
      for (i = 0; i < imageSize[box]; i++) {
        cosMolRef[*thisMol][i] = 0.0;
        sinMolRef[*thisMol][i] = 0.0;

        for (j = 0; j < thisKind.NumAtoms(); j++) {
          dotProduct = currentAxes.DotProduct(mols.MolStart(*thisMol) + j,
                                              kx[box][i], ky[box][i],
                                              kz[box][i], molCoords);

          cosMolRef[*thisMol][i] += (thisKind.AtomCharge(j) *
                                     cos(dotProduct));
          sinMolRef[*thisMol][i] += (thisKind.AtomCharge(j) *
                                     sin(dotProduct));
        }
        sumRnew[box][i] += cosMolRef[*thisMol][i];
        sumInew[box][i] += sinMolRef[*thisMol][i];
      }
      thisMol++;
    }
  }
}


//calculate reciprocate term for a box
double EwaldCached::BoxReciprocal(uint box) const
{
  int i;
  double energyRecip = 0.0;

  if (box < BOXES_WITH_U_NB) {
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i) reduction(+:energyRecip)
#endif
    for (i = 0; i < imageSize[box]; i++) {
      energyRecip += (( sumRnew[box][i] * sumRnew[box][i] +
                        sumInew[box][i] * sumInew[box][i]) *
                      prefact[box][i]);
    }
  }

  return energyRecip;
}

//calculate reciprocate term for displacement and rotation move
double EwaldCached::MolReciprocal(XYZArray const& molCoords,
                                  const uint molIndex,
                                  const uint box,
                                  XYZ const*const newCOM)
{
  double energyRecipNew = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const& thisKind = mols.GetKind(molIndex);
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);
    uint p, atom;
	int i;
    double sumRealNew, sumImaginaryNew, dotProductNew, sumRealOld,
           sumImaginaryOld;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, atom, sumRealNew, sumImaginaryNew, sumRealOld, sumImaginaryOld, dotProductNew) reduction(+:energyRecipNew)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      sumRealNew = 0.0;
      sumImaginaryNew = 0.0;
      dotProductNew = 0.0;
      sumRealOld = cosMolRef[molIndex][i];
      sumImaginaryOld = sinMolRef[molIndex][i];
      cosMolRestore[i] = cosMolRef[molIndex][i];
      sinMolRestore[i] = sinMolRef[molIndex][i];

      for (p = 0; p < length; ++p) {
        atom = startAtom + p;
        dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
                                               kyRef[box][i], kzRef[box][i],
                                               molCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));
      }

      sumRnew[box][i] = sumRref[box][i] - sumRealOld + sumRealNew;
      sumInew[box][i] = sumIref[box][i] - sumImaginaryOld + sumImaginaryNew;
      cosMolRef[molIndex][i] = sumRealNew;
      sinMolRef[molIndex][i] = sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
  }

  return energyRecipNew - sysPotRef.boxEnergy[box].recip;
}

//calculate reciprocate term in destination box for swap move
double EwaldCached::SwapDestRecip(const cbmc::TrialMol &newMol,
                                  const uint box, const int sourceBox,
                                  const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

#ifdef _OPENMP
  #pragma omp parallel default(shared)
#endif
  {
    std::memcpy(cosMolRestore, cosMolRef[molIndex], sizeof(double)*imageLarge);
    std::memcpy(sinMolRestore, sinMolRef[molIndex], sizeof(double)*imageLarge);
  }

  if (box < BOXES_WITH_U_NB) {
    uint p, length;
	int i;
    MoleculeKind const& thisKind = newMol.GetKind();
    XYZArray molCoords = newMol.GetCoords();
    double dotProductNew;
    length = thisKind.NumAtoms();

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, dotProductNew) reduction(+:energyRecipNew)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      cosMolRef[molIndex][i] = 0.0;
      sinMolRef[molIndex][i] = 0.0;
      dotProductNew = 0.0;

      for (p = 0; p < length; ++p) {
        dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
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


//calculate reciprocate term in source box for swap move
double EwaldCached::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                    const uint box, const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

  if (box < BOXES_WITH_U_NB) {
    int i;
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i) reduction(+:energyRecipNew)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      sumRnew[box][i] = sumRref[box][i] - cosMolRestore[i];
      sumInew[box][i] = sumIref[box][i] - sinMolRestore[i];

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }

    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }
  return energyRecipNew - energyRecipOld;
}

void EwaldCached::UpdateRecipVec(uint box)
{
  double *tempKx, *tempKy, *tempKz, *tempHsqr, *tempPrefact;
  tempKx = kxRef[box];
  tempKy = kyRef[box];
  tempKz = kzRef[box];
  tempHsqr = hsqrRef[box];
  tempPrefact = prefactRef[box];

  kxRef[box] = kx[box];
  kyRef[box] = ky[box];
  kzRef[box] = kz[box];
  hsqrRef[box] = hsqr[box];
  prefactRef[box] = prefact[box];

  kx[box] = tempKx;
  ky[box] = tempKy;
  kz[box] = tempKz;
  hsqr[box] = tempHsqr;
  prefact[box] = tempPrefact;
#ifdef GOMC_CUDA
  UpdateRecipVecCUDA(forcefield.particles->getCUDAVars(), box);
#endif

  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    imageSizeRef[b] = imageSize[b];
  }
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
