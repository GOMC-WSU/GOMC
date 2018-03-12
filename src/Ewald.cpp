/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Ewald.h"
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
#include "CalculateEwaldCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
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

Ewald::Ewald(StaticVals & stat, System & sys) :
  forcefield(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
  currentCOM(sys.com), sysPotRef(sys.potential),
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup(sys.molLookup),
#else
  molLookup(stat.molLookup),
#endif
#ifdef VARIABLE_VOLUME
  currentAxes(sys.boxDimRef)
#else
  currentAxes(*stat.GetBoxDim())
#endif
{
  ewald = false;
  electrostatic = false;
  imageLarge = 0;
  alpha = 0.0;
  recip_rcut = 0.0;
  recip_rcut_Sq = 0.0;
}

Ewald::~Ewald()
{
  if(ewald) {
#ifdef GOMC_CUDA
    DestroyEwaldCUDAVars(forcefield.particles->getCUDAVars());
#endif
    for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
      if (kx[b] != NULL) {
        delete[] kx[b];
        delete[] ky[b];
        delete[] kz[b];
        delete[] hsqr[b];
        delete[] prefact[b];
        delete[] kxRef[b];
        delete[] kyRef[b];
        delete[] kzRef[b];
        delete[] hsqrRef[b];
        delete[] prefactRef[b];
        delete[] sumRnew[b];
        delete[] sumInew[b];
        delete[] sumRref[b];
        delete[] sumIref[b];
      }
    }

    if (kx != NULL) {
      delete[] kmax;
      delete[] kx;
      delete[] ky;
      delete[] kz;
      delete[] hsqr;
      delete[] prefact;
      delete[] kxRef;
      delete[] kyRef;
      delete[] kzRef;
      delete[] hsqrRef;
      delete[] prefactRef;
      delete[] sumRnew;
      delete[] sumInew;
      delete[] sumRref;
      delete[] sumIref;
      delete[] imageSize;
      delete[] imageSizeRef;
    }
  }
}

void Ewald::Init()
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
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    RecipInit(b, currentAxes);
    BoxReciprocalSetup(b, currentCoords);
    SetRecipRef(b);
    printf("Box: %d, RecipVectors: %d, kmax: %d\n", b, imageSize[b],
           kmax[b]);
  }
}

void Ewald::AllocMem()
{
  //get size of image using defined Kmax
  //Allocate Memory

  kmax = new uint[BOXES_WITH_U_NB];
  imageSize = new uint[BOXES_WITH_U_NB];
  imageSizeRef = new uint[BOXES_WITH_U_NB];
  sumRnew = new double*[BOXES_WITH_U_NB];
  sumInew = new double*[BOXES_WITH_U_NB];
  sumRref = new double*[BOXES_WITH_U_NB];
  sumIref = new double*[BOXES_WITH_U_NB];
  kx = new double*[BOXES_WITH_U_NB];
  ky = new double*[BOXES_WITH_U_NB];
  kz = new double*[BOXES_WITH_U_NB];
  hsqr = new double*[BOXES_WITH_U_NB];
  prefact = new double*[BOXES_WITH_U_NB];
  kxRef = new double*[BOXES_WITH_U_NB];
  kyRef = new double*[BOXES_WITH_U_NB];
  kzRef = new double*[BOXES_WITH_U_NB];
  hsqrRef = new double*[BOXES_WITH_U_NB];
  prefactRef = new double*[BOXES_WITH_U_NB];

  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    RecipCountInit(b, currentAxes);
  }
  //25% larger than original box size, reserved for image size change
  imageTotal = findLargeImage();

  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    kx[b] = new double[imageTotal];
    ky[b] = new double[imageTotal];
    kz[b] = new double[imageTotal];
    hsqr[b] = new double[imageTotal];
    prefact[b] = new double[imageTotal];
    kxRef[b] = new double[imageTotal];
    kyRef[b] = new double[imageTotal];
    kzRef[b] = new double[imageTotal];
    hsqrRef[b] = new double[imageTotal];
    prefactRef[b] = new double[imageTotal];
    sumRnew[b] = new double[imageTotal];
    sumInew[b] = new double[imageTotal];
    sumRref[b] = new double[imageTotal];
    sumIref[b] = new double[imageTotal];
  }

#ifdef GOMC_CUDA
  InitEwaldVariablesCUDA(forcefield.particles->getCUDAVars(), imageTotal);
#endif
}


//calculate reciprocate term for a box
void Ewald::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
  uint j, m;
  int i;
  double dotProduct = 0.0;
  double sumReal = 0.0;
  double sumImaginary = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);

#ifdef GOMC_CUDA
    int numberOfAtoms = 0, i = 0;

    for(int k = 0; k < mols.GetKindsCount(); k++) {
      MoleculeKind const& thisKind = mols.kinds[k];
      numberOfAtoms += thisKind.NumAtoms() * molLookup.NumKindInBox(k, box);
    }

    XYZArray thisBoxCoords(numberOfAtoms);
    std::vector<double> chargeBox;

    while (thisMol != end) {
      MoleculeKind const& thisKind = mols.GetKind(*thisMol);
      for (j = 0; j < thisKind.NumAtoms(); j++) {
        thisBoxCoords.Set(i, molCoords[mols.MolStart(*thisMol) + j]);
        chargeBox.push_back(thisKind.AtomCharge(j));
        i++;
      }
      thisMol++;
    }
    CallBoxReciprocalSetupGPU(forcefield.particles->getCUDAVars(),
                              thisBoxCoords, kx[box], ky[box], kz[box],
                              chargeBox, imageSize[box], sumRnew[box],
                              sumInew[box], prefact[box], hsqr[box],
                              currentEnergyRecip[box], box);
#else
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
      #pragma omp parallel for default(shared) private(i, j, dotProduct, sumReal, sumImaginary)
#endif
      for (i = 0; i < imageSize[box]; i++) {
        sumReal = 0.0;
        sumImaginary = 0.0;

        for (j = 0; j < thisKind.NumAtoms(); j++) {
          dotProduct = currentAxes.DotProduct(mols.MolStart(*thisMol) + j,
                                              kx[box][i], ky[box][i],
                                              kz[box][i], molCoords);

          sumReal += (thisKind.AtomCharge(j) * cos(dotProduct));
          sumImaginary += (thisKind.AtomCharge(j) * sin(dotProduct));
        }
        sumRnew[box][i] += sumReal;
        sumInew[box][i] += sumImaginary;
      }
      thisMol++;
    }
#endif
  }
}


//calculate reciprocate term for a box
double Ewald::BoxReciprocal(uint box) const
{
  int i;
  double energyRecip = 0.0;

  if (box < BOXES_WITH_U_NB) {
#ifdef GOMC_CUDA
    return currentEnergyRecip[box];
#else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i) reduction(+:energyRecip)
#endif
    for (i = 0; i < imageSize[box]; i++) {
      energyRecip += (( sumRnew[box][i] * sumRnew[box][i] +
                        sumInew[box][i] * sumInew[box][i]) *
                      prefact[box][i]);
    }
#endif
  }

  return energyRecip;
}


//calculate reciprocate term for displacement and rotation move
double Ewald::MolReciprocal(XYZArray const& molCoords,
                            const uint molIndex, const uint box,
                            XYZ const*const newCOM)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const& thisKind = mols.GetKind(molIndex);
    uint length = thisKind.NumAtoms();
    uint startAtom = mols.MolStart(molIndex);
    uint p, atom;
	int i;
    double sumRealNew, sumImaginaryNew, dotProductNew, dotProductOld,
           sumRealOld, sumImaginaryOld;
#ifdef GOMC_CUDA
    XYZArray cCoords(length);
    std::vector<double> MolCharge;
    for(p = 0; p < length; p++) {
      cCoords.Set(p, currentCoords[startAtom + p]);
      MolCharge.push_back(thisKind.AtomCharge(p));
    }
    CallMolReciprocalGPU(forcefield.particles->getCUDAVars(),
                         cCoords, molCoords, MolCharge, imageSizeRef[box],
                         sumRnew[box], sumInew[box], energyRecipNew, box);
#else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, atom, sumRealNew, sumImaginaryNew, sumRealOld, sumImaginaryOld, dotProductNew, dotProductOld) reduction(+:energyRecipNew, energyRecipOld)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      sumRealNew = 0.0;
      sumImaginaryNew = 0.0;
      dotProductNew = 0.0;
      dotProductOld = 0.0;
      sumRealOld = 0.0;
      sumImaginaryOld = 0.0;

      for (p = 0; p < length; ++p) {
        atom = startAtom + p;
        dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
                                               kyRef[box][i], kzRef[box][i],
                                               molCoords);

        dotProductOld = currentAxes.DotProduct(atom, kxRef[box][i],
                                               kyRef[box][i], kzRef[box][i],
                                               currentCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));

        sumRealOld += (thisKind.AtomCharge(p) * cos(dotProductOld));
        sumImaginaryOld += (thisKind.AtomCharge(p) * sin(dotProductOld));
      }

      sumRnew[box][i] = sumRref[box][i] - sumRealOld + sumRealNew;
      sumInew[box][i] = sumIref[box][i] - sumImaginaryOld + sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
#endif
  }
  return energyRecipNew - sysPotRef.boxEnergy[box].recip;
}



//calculate reciprocate term in destination box for swap move
double Ewald::SwapDestRecip(const cbmc::TrialMol &newMol,
                            const uint box, const int sourceBox,
                            const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;


  if (box < BOXES_WITH_U_NB) {
    uint p, length;
	int i;
    MoleculeKind const& thisKind = newMol.GetKind();
    XYZArray molCoords = newMol.GetCoords();
    double dotProductNew, sumRealNew, sumImaginaryNew;
    length = thisKind.NumAtoms();
#ifdef GOMC_CUDA
    bool insert = true;
    std::vector<double> MolCharge;
    for(p = 0; p < length; p++) {
      MolCharge.push_back(thisKind.AtomCharge(p));
    }
    CallSwapReciprocalGPU(forcefield.particles->getCUDAVars(),
                          molCoords, MolCharge, imageSizeRef[box],
                          sumRnew[box], sumInew[box],
                          insert, energyRecipNew, box);
#else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, dotProductNew, sumRealNew, sumImaginaryNew) reduction(+:energyRecipNew)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      sumRealNew = 0.0;
      sumImaginaryNew = 0.0;
      dotProductNew = 0.0;

      for (p = 0; p < length; ++p) {
        dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
                                               kyRef[box][i], kzRef[box][i],
                                               molCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));
      }

      //sumRealNew;
      sumRnew[box][i] = sumRref[box][i] + sumRealNew;
      //sumImaginaryNew;
      sumInew[box][i] = sumIref[box][i] + sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
#endif
    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }

  return energyRecipNew - energyRecipOld;
}

void Ewald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  if(boxAxes.orthogonal[box])
    RecipInitOrth(box, boxAxes);
  else
    RecipInitNonOrth(box, boxAxes);
}


//calculate reciprocate term in source box for swap move
double Ewald::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                              const uint box, const int molIndex)
{
  double energyRecipNew = 0.0;
  double energyRecipOld = 0.0;

  if (box < BOXES_WITH_U_NB) {
    uint p;
	int i;
    double sumRealNew, sumImaginaryNew, dotProductNew;
    MoleculeKind const& thisKind = oldMol.GetKind();
    XYZArray molCoords = oldMol.GetCoords();
    uint length = thisKind.NumAtoms();
#ifdef GOMC_CUDA
    bool insert = false;
    std::vector<double> MolCharge;
    for(p = 0; p < length; p++) {
      MolCharge.push_back(thisKind.AtomCharge(p));
    }
    CallSwapReciprocalGPU(forcefield.particles->getCUDAVars(),
                          molCoords, MolCharge, imageSizeRef[box],
                          sumRnew[box], sumInew[box],
                          insert, energyRecipNew, box);

#else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, dotProductNew, sumRealNew, sumImaginaryNew) reduction(+:energyRecipNew)
#endif
    for (i = 0; i < imageSizeRef[box]; i++) {
      sumRealNew = 0.0;
      sumImaginaryNew = 0.0;
      dotProductNew = 0.0;

      for (p = 0; p < length; ++p) {
        dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
                                               kyRef[box][i], kzRef[box][i],
                                               molCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));

      }
      sumRnew[box][i] = sumRref[box][i] - sumRealNew;
      sumInew[box][i] = sumIref[box][i] - sumImaginaryNew;

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
#endif
    energyRecipOld = sysPotRef.boxEnergy[box].recip;
  }
  return energyRecipNew - energyRecipOld;
}


//restore cosMol and sinMol
void Ewald::RestoreMol(int molIndex)
{
  return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Ewald::exgMolCache()
{
  return;
}

//compare number of images in different boxes and select the largest one
uint Ewald::findLargeImage()
{
  imageLarge = 0;
  for (int b = 0; b < BOXES_WITH_U_NB; b++) {
    if (imageLarge < imageSize[b])
      imageLarge = imageSize[b];
  }
  return imageLarge;
}

void Ewald::RecipInitOrth(uint box, BoxDimensions const& boxAxes)
{
  uint counter = 0;
  int x, y, z, nkx_max, nky_max, nky_min, nkz_max, nkz_min;
  double ksqr, kX, kY, kZ;
  double alpsqr4 = 1.0 / (4.0 * alpha * alpha);
  XYZ constValue = boxAxes.axis.Get(box);
  constValue.Inverse();
  constValue *= 2 * M_PI;

  double vol = boxAxes.volume[box] / (4 * M_PI);
  nkx_max = int(recip_rcut * boxAxes.axis.Get(box).x / (2 * M_PI)) + 1;
  nky_max = int(recip_rcut * boxAxes.axis.Get(box).y / (2 * M_PI)) + 1;
  nkz_max = int(recip_rcut * boxAxes.axis.Get(box).z / (2 * M_PI)) + 1;
  kmax[box] = std::max(std::max(nkx_max, nky_max), std::max(nky_max, nkz_max));

  for(x = 0; x <= nkx_max; x++) {
    if(x == 0.0)
      nky_min = 0;
    else
      nky_min = -nky_max;

    for(y = nky_min; y <= nky_max; y++) {
      if(x == 0.0 && y == 0.0)
        nkz_min = 1;
      else
        nkz_min = -nkz_max;

      for(z = nkz_min; z <= nkz_max; z++) {
        kX = constValue.x * x;
        kY = constValue.y * y;
        kZ = constValue.z * z;
        ksqr = kX * kX + kY * kY + kZ * kZ;

        if(ksqr < recip_rcut_Sq) {
          kx[box][counter] = kX;
          ky[box][counter] = kY;
          kz[box][counter] = kZ;
          hsqr[box][counter] = ksqr;
          prefact[box][counter] = num::qqFact * exp(-ksqr * alpsqr4) /
                                  (ksqr * vol);
          counter++;
        }
      }
    }
  }

  imageSize[box] = counter;

  if (counter > imageTotal) {
    std::cout << "Error: Kmax exceeded due to large change in system volume.\n";
    std::cout << "Restart the simulation from restart files.\n";
    exit(EXIT_FAILURE);
  }
}

void Ewald::RecipInitNonOrth(uint box, BoxDimensions const& boxAxes)
{
  uint counter = 0;
  int x, y, z, nkx_max, nky_max, nky_min, nkz_max, nkz_min;
  double ksqr, kX, kY, kZ;
  double alpsqr4 = 1.0 / (4.0 * alpha * alpha);
  XYZArray cellB(boxAxes.cellBasis[box]);
  cellB.Scale(0, boxAxes.axis.Get(box).x);
  cellB.Scale(1, boxAxes.axis.Get(box).y);
  cellB.Scale(2, boxAxes.axis.Get(box).z);
  XYZArray cellB_Inv(3);
  double det = cellB.AdjointMatrix(cellB_Inv);
  cellB_Inv.ScaleRange(0, 3, (2 * M_PI) / det);

  double vol = boxAxes.volume[box] / (4 * M_PI);
  nkx_max = int(recip_rcut * boxAxes.axis.Get(box).x / (2 * M_PI)) + 1;
  nky_max = int(recip_rcut * boxAxes.axis.Get(box).y / (2 * M_PI)) + 1;
  nkz_max = int(recip_rcut * boxAxes.axis.Get(box).z / (2 * M_PI)) + 1;
  kmax[box] = std::max(std::max(nkx_max, nky_max), std::max(nky_max, nkz_max));

  for (x = 0; x <= nkx_max; x++) {
    if(x == 0.0)
      nky_min = 0;
    else
      nky_min = -nky_max;

    for(y = nky_min; y <= nky_max; y++) {
      if(x == 0.0 && y == 0.0)
        nkz_min = 1;
      else
        nkz_min = -nkz_max;

      for(z = nkz_min; z <= nkz_max; z++) {
        kX = boxAxes.DotProduct(cellB_Inv.Get(0), XYZ(x, y, z));
        kY = boxAxes.DotProduct(cellB_Inv.Get(1), XYZ(x, y, z));
        kZ = boxAxes.DotProduct(cellB_Inv.Get(2), XYZ(x, y, z));
        ksqr = kX * kX + kY * kY + kZ * kZ;

        if(ksqr < recip_rcut_Sq) {
          kx[box][counter] = kX;
          ky[box][counter] = kY;
          kz[box][counter] = kZ;
          hsqr[box][counter] = ksqr;
          prefact[box][counter] = num::qqFact * exp(-ksqr * alpsqr4) /
                                  (ksqr * vol);
          counter++;
        }
      }
    }
  }

  imageSize[box] = counter;

  if(counter > imageTotal) {
    std::cout << "Error: Kmax exceeded due to large change in system volume.\n";
    std::cout << "Restart the simulation from restart files.\n";
    exit(EXIT_FAILURE);
  }
}

//estimate number of vectors
void Ewald::RecipCountInit(uint box, BoxDimensions const& boxAxes)
{
  uint counter = 0;
  int x, y, z, nkx_max, nky_max, nky_min, nkz_max, nkz_min;
  double ksqr, excess, kX, kY, kZ;
  XYZArray cellB(boxAxes.cellBasis[box]);
  XYZ constValue = boxAxes.axis.Get(box);
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  excess = 1.25;
#else
  excess = 1.00;
#endif
  constValue *= excess;
  cellB.Scale(0, constValue.x);
  cellB.Scale(1, constValue.y);
  cellB.Scale(2, constValue.z);
  XYZArray cellB_Inv(3);
  double det = cellB.AdjointMatrix(cellB_Inv);
  cellB_Inv.ScaleRange(0, 3, (2 * M_PI) / det);

  nkx_max = int(recip_rcut * constValue.x / (2 * M_PI)) + 1;
  nky_max = int(recip_rcut * constValue.y / (2 * M_PI)) + 1;
  nkz_max = int(recip_rcut * constValue.z / (2 * M_PI)) + 1;

  for(x = 0; x <= nkx_max; x++) {
    if(x == 0.0)
      nky_min = 0;
    else
      nky_min = -nky_max;

    for(y = nky_min; y <= nky_max; y++) {
      if(x == 0.0 && y == 0.0)
        nkz_min = 1;
      else
        nkz_min = -nkz_max;

      for(z = nkz_min; z <= nkz_max; z++) {
        kX = boxAxes.DotProduct(cellB_Inv.Get(0), XYZ(x, y, z));
        kY = boxAxes.DotProduct(cellB_Inv.Get(1), XYZ(x, y, z));
        kZ = boxAxes.DotProduct(cellB_Inv.Get(2), XYZ(x, y, z));
        ksqr = kX * kX + kY * kY + kZ * kZ;

        if(ksqr < recip_rcut_Sq) {
          counter++;
        }
      }
    }
  }
  imageSize[box] = counter;
}

//back up reciptocate value to Ref (will be called during initialization)
void Ewald::SetRecipRef(uint box)
{
#ifdef _OPENMP
  #pragma omp parallel default(shared)
#endif
  {
    std::memcpy(sumRref[box], sumRnew[box], sizeof(double) * imageSize[box]);
    std::memcpy(sumIref[box], sumInew[box], sizeof(double) * imageSize[box]);
    std::memcpy(kxRef[box], kx[box], sizeof(double) * imageSize[box]);
    std::memcpy(kyRef[box], ky[box], sizeof(double) * imageSize[box]);
    std::memcpy(kzRef[box], kz[box], sizeof(double) * imageSize[box]);
    std::memcpy(hsqrRef[box], hsqr[box], sizeof(double) * imageSize[box]);
    std::memcpy(prefactRef[box], prefact[box], sizeof(double) *imageSize[box]);
  }
#ifdef GOMC_CUDA
  CopyCurrentToRefCUDA(forcefield.particles->getCUDAVars(),
                       box, imageSize[box]);
#endif
  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    imageSizeRef[b] = imageSize[b];
  }
}