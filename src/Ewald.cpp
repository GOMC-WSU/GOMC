/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Ewald.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"             //
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "GeomLib.h"
#include "NumLib.h"
#include "pme.h"
#include <cassert>
#ifdef GOMC_CUDA
#include "CalculateEwaldCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
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
  ff(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
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
  alpha = 0.0;
  recip_rcut = 0.0;
  recip_rcut_Sq = 0.0;
  multiParticleEnabled = stat.multiParticleEnabled;
}

Ewald::~Ewald()
{
  if(ff.ewald) {
#ifdef GOMC_CUDA
    DestroyEwaldCUDAVars(ff.particles->getCUDAVars());
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
      if(abs(molKind.AtomCharge(a)) < 0.000000001) {
        particleHasNoCharge.push_back(true);
      } else {
        particleHasNoCharge.push_back(false);
      }
    }
  }

  AllocMem();
  //initialize K vectors and reciprocate terms
  UpdateVectorsAndRecipTerms();
}

void Ewald::UpdateVectorsAndRecipTerms()
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    RecipInit(b, currentAxes);
    BoxReciprocalSetup(b, currentCoords);
    SetRecipRef(b);
    printf("Box: %d, RecipVectors: %6d, kmax: %d\n",
           b, imageSize[b], kmax[b]);
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
  InitEwaldVariablesCUDA(ff.particles->getCUDAVars(), imageTotal);
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
    CallBoxReciprocalSetupGPU(ff.particles->getCUDAVars(),
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
      uint start = mols.MolStart(*thisMol);

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(i, j, dotProduct, sumReal, sumImaginary)
#endif
      for (i = 0; i < imageSize[box]; i++) {
        sumReal = 0.0;
        sumImaginary = 0.0;

        for (j = 0; j < thisKind.NumAtoms(); j++) {
          if(particleHasNoCharge[start + j]) {
            continue;
          }
          dotProduct = Dot(start + j, kx[box][i], ky[box][i],
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
                            const uint molIndex, const uint box)
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
    CallMolReciprocalGPU(ff.particles->getCUDAVars(),
                         cCoords, molCoords, MolCharge, imageSizeRef[box],
                         sumRnew[box], sumInew[box], energyRecipNew, box);
#else
#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, atom, sumRealNew, \
    sumImaginaryNew, sumRealOld, sumImaginaryOld, dotProductNew, dotProductOld) \
reduction(+:energyRecipNew, energyRecipOld)
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
        if(particleHasNoCharge[atom]) {
          continue;
        }
        dotProductNew = Dot(p, kxRef[box][i],
                            kyRef[box][i], kzRef[box][i],
                            molCoords);

        dotProductOld = Dot(atom, kxRef[box][i],
                            kyRef[box][i], kzRef[box][i],
                            currentCoords);

        sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
        sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));

        sumRealOld += (thisKind.AtomCharge(p) * cos(dotProductOld));
        sumImaginaryOld += (thisKind.AtomCharge(p) * sin(dotProductOld));
      }

      sumRnew[box][i] = sumRref[box][i] + (sumRealNew - sumRealOld);
      sumInew[box][i] = sumIref[box][i] + (sumImaginaryNew - sumImaginaryOld);

      energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
                         * sumInew[box][i]) * prefactRef[box][i];
    }
#endif
  }
  return energyRecipNew - sysPotRef.boxEnergy[box].recip;
}

void Ewald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  if(boxAxes.orthogonal[box])
    RecipInitOrth(box, boxAxes);
  else
    RecipInitNonOrth(box, boxAxes);
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
  uint maxImg = 0;
  for (int b = 0; b < BOXES_WITH_U_NB; b++) {
    if (maxImg < imageSize[b])
      maxImg = imageSize[b];
  }
  return maxImg;
}

//backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Ewald::backupMolCache()
{
  return;
}

void Ewald::RecipInitOrth(uint box, BoxDimensions const& boxAxes)
{
  uint counter = 0;
  int x, y, z, nkx_max, nky_max, nky_min, nkz_max, nkz_min;
  double ksqr, kX, kY, kZ;
  double alpsqr4 = 1.0 / (4.0 * ff.alphaSq[box]);
  XYZ constValue = boxAxes.axis.Get(box);
  constValue.Inverse();
  constValue *= 2 * M_PI;

  double vol = boxAxes.volume[box] / (4 * M_PI);
  nkx_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).x / (2 * M_PI)) + 1;
  nky_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).y / (2 * M_PI)) + 1;
  nkz_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).z / (2 * M_PI)) + 1;
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

        if(ksqr < ff.recip_rcut_Sq[box]) {
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
  double alpsqr4 = 1.0 / (4.0 * ff.alphaSq[box]);
  XYZArray cellB(boxAxes.cellBasis[box]);
  cellB.Scale(0, boxAxes.axis.Get(box).x);
  cellB.Scale(1, boxAxes.axis.Get(box).y);
  cellB.Scale(2, boxAxes.axis.Get(box).z);
  XYZArray cellB_Inv(3);
  double det = cellB.AdjointMatrix(cellB_Inv);
  cellB_Inv.ScaleRange(0, 3, (2 * M_PI) / det);

  double vol = boxAxes.volume[box] / (4 * M_PI);
  nkx_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).x / (2 * M_PI)) + 1;
  nky_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).y / (2 * M_PI)) + 1;
  nkz_max = int(ff.recip_rcut[box] * boxAxes.axis.Get(box).z / (2 * M_PI)) + 1;
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
        kX = Dot(cellB_Inv.Get(0), XYZ(x, y, z));
        kY = Dot(cellB_Inv.Get(1), XYZ(x, y, z));
        kZ = Dot(cellB_Inv.Get(2), XYZ(x, y, z));
        ksqr = kX * kX + kY * kY + kZ * kZ;

        if(ksqr < ff.recip_rcut_Sq[box]) {
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
#if ENSEMBLE == GEMC
  excess = 1.25;
#elif ENSEMBLE == NPT
  excess = 1.5;
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

  nkx_max = int(ff.recip_rcut[box] * constValue.x / (2 * M_PI)) + 1;
  nky_max = int(ff.recip_rcut[box] * constValue.y / (2 * M_PI)) + 1;
  nkz_max = int(ff.recip_rcut[box] * constValue.z / (2 * M_PI)) + 1;

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
        kX = Dot(cellB_Inv.Get(0), XYZ(x, y, z));
        kY = Dot(cellB_Inv.Get(1), XYZ(x, y, z));
        kZ = Dot(cellB_Inv.Get(2), XYZ(x, y, z));
        ksqr = kX * kX + kY * kY + kZ * kZ;

        if(ksqr < ff.recip_rcut_Sq[box]) {
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
  CopyCurrentToRefCUDA(ff.particles->getCUDAVars(),
                       box, imageSize[box]);
#endif
  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    imageSizeRef[b] = imageSize[b];
  }
}

double Ewald::MolCorrection(uint molIndex, uint box) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  double dist, distSq;
  double correction = 0.0;
  XYZ virComponents;

  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);

  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box);
      dist = sqrt(distSq);
      correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                     erf(ff.alpha[box] * dist) / dist);
    }
  }

  return correction;
}

//calculate self term for a box
double Ewald::BoxSelf(BoxDimensions const& boxAxes, uint box) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  double self = 0.0;
  double molSelfEnergy;
  uint i, j, length;

  for (i = 0; i < mols.GetKindsCount(); i++) {
    MoleculeKind const& thisKind = mols.kinds[i];
    length = thisKind.NumAtoms();
    molSelfEnergy = 0.0;

    for (j = 0; j < length; j++) {
      molSelfEnergy += (thisKind.AtomCharge(j) * thisKind.AtomCharge(j));
    }
    self += (molSelfEnergy * molLookup.NumKindInBox(i, box));
  }

  self = -1.0 * self * ff.alpha[box] * num::qqFact / sqrt(M_PI);

  return self;
}

// NOTE: The calculation of W12, W13, W23 is expensive and would not be
// requied for pressure and surface tension calculation. So, they have been
// commented out. In case you need to calculate them, uncomment them.
Virial Ewald::VirialReciprocal(Virial& virial, uint box) const
{
  Virial tempVir = virial;
  if (box >= BOXES_WITH_U_NB)
    return tempVir;

  double wT11 = 0.0, wT12 = 0.0, wT13 = 0.0;
  double wT22 = 0.0, wT23 = 0.0, wT33 = 0.0;

  double recipIntra = 0.0;
  double constVal = 1.0 / (4.0 * ff.alphaSq[box]);
  double factor, arg, charge;
  uint p, length, start, atom;
  int i;

  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box),
                               end = molLookup.BoxEnd(box);

  XYZ atomC, comC, diffC;

#ifdef GOMC_CUDA
  int numberOfAtoms = 0, atomIndex = 0;

  for(int k = 0; k < mols.GetKindsCount(); k++) {
    MoleculeKind const& thisKind = mols.kinds[k];
    numberOfAtoms += thisKind.NumAtoms() * molLookup.NumKindInBox(k, box);
  }

  XYZArray thisBoxCoords(numberOfAtoms);
  XYZArray thisBoxCOMDiff(numberOfAtoms);
  std::vector<double> chargeBox;

  while (thisMol != end) {
    length = mols.GetKind(*thisMol).NumAtoms();
    start = mols.MolStart(*thisMol);
    comC = currentCOM.Get(*thisMol);

    for (p = 0; p < length; p++) {
      atom = start + p;
      //compute the vector of the bead to the COM (p)
      // need to unwrap the atom coordinate
      atomC = currentCoords.Get(atom);
      currentAxes.UnwrapPBC(atomC, box, comC);
      diffC = atomC - comC;

      thisBoxCoords.Set(atomIndex, atomC);
      thisBoxCOMDiff.Set(atomIndex, diffC);
      chargeBox.push_back(particleCharge[atom]);
      atomIndex++;
    }
    thisMol++;
  }

  CallVirialReciprocalGPU(ff.particles->getCUDAVars(), thisBoxCoords,
                          thisBoxCOMDiff, chargeBox, wT11, wT12,
                          wT13, wT22, wT23, wT33, imageSizeRef[box], constVal,
                          box);
#else
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(i, factor) reduction(+:wT11, wT22, wT33)
#endif
  for (i = 0; i < imageSizeRef[box]; i++) {
    factor = prefactRef[box][i] * (sumRref[box][i] * sumRref[box][i] +
                                   sumIref[box][i] * sumIref[box][i]);

    wT11 += factor * (1.0 - 2.0 * (constVal + 1.0 / hsqrRef[box][i]) *
                      kxRef[box][i] * kxRef[box][i]);

    wT22 += factor * (1.0 - 2.0 * (constVal + 1.0 / hsqrRef[box][i]) *
                      kyRef[box][i] * kyRef[box][i]);

    wT33 += factor * (1.0 - 2.0 * (constVal + 1.0 / hsqrRef[box][i]) *
                      kzRef[box][i] * kzRef[box][i]);
  }

  //Intramolecular part
  while (thisMol != end) {
    length = mols.GetKind(*thisMol).NumAtoms();
    start = mols.MolStart(*thisMol);
    comC = currentCOM.Get(*thisMol);

    for (p = 0; p < length; p++) {
      atom = start + p;
      if(particleHasNoCharge[atom]) {
        continue;
      }
      //compute the vector of the bead to the COM (p)
      // need to unwrap the atom coordinate
      atomC = currentCoords.Get(atom);
      currentAxes.UnwrapPBC(atomC, box, comC);
      diffC = atomC - comC;
      charge = particleCharge[atom];

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(i, arg, factor) reduction(+:wT11, wT22, wT33)
#endif
      for (i = 0; i < imageSizeRef[box]; i++) {
        //compute the dot product of k and r
        arg = Dot(atom, kxRef[box][i], kyRef[box][i],
                  kzRef[box][i], currentCoords);

        factor = prefactRef[box][i] * 2.0 * (sumIref[box][i] * cos(arg) -
                                             sumRref[box][i] * sin(arg)) *
                 charge;

        wT11 += factor * (kxRef[box][i] * diffC.x);

        wT22 += factor * (kyRef[box][i] * diffC.y);

        wT33 += factor * (kzRef[box][i] * diffC.z);
      }
    }
    ++thisMol;
  }
#endif

  // set the all tensor values
  tempVir.recipTens[0][0] = wT11;
  tempVir.recipTens[0][1] = wT12;
  tempVir.recipTens[0][2] = wT13;

  tempVir.recipTens[1][0] = wT12;
  tempVir.recipTens[1][1] = wT22;
  tempVir.recipTens[1][2] = wT23;

  tempVir.recipTens[2][0] = wT13;
  tempVir.recipTens[2][1] = wT23;
  tempVir.recipTens[2][2] = wT33;

  // setting virial of reciprocal cpace
  tempVir.recip = wT11 + wT22 + wT33;

  return tempVir;
}

//update reciprocate values
void Ewald::UpdateRecip(uint box)
{
  double *tempR, *tempI;
  tempR = sumRref[box];
  tempI = sumIref[box];
  sumRref[box] = sumRnew[box];
  sumIref[box] = sumInew[box];
  sumRnew[box] = tempR;
  sumInew[box] = tempI;
#ifdef GOMC_CUDA
  UpdateRecipCUDA(ff.particles->getCUDAVars(), box);
#endif
}

void Ewald::UpdateRecipVec(uint box)
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
  UpdateRecipVecCUDA(ff.particles->getCUDAVars(), box);
#endif

  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    imageSizeRef[b] = imageSize[b];
  }
}


//calculate reciprocate force term for a box with molCoords
double Ewald::BoxForceReciprocal(XYZArray const& molCoords,
                               XYZArray& atomForceRec,
                               XYZArray& molForceRec,
                               uint box)
{
  if(multiParticleEnabled && (box < BOXES_WITH_U_NB)) {
    double LL[3 * 3];
    XYZ boxDimension = currentAxes.axis.Get(box);
    for(int i=0; i<3; i++) {
      LL[3*i + 0] = currentAxes.cellBasis[box][i].x;
      LL[3*i + 1] = currentAxes.cellBasis[box][i].y;
      LL[3*i + 2] = currentAxes.cellBasis[box][i].z;
    }

    cout << "currentAxes.rCut: " << currentAxes.rCut[box] << endl;
    ewald::pme p(boxDimension.x, boxDimension.y, boxDimension.z, LL, 
                 currentCoords.Count(), particleCharge.data(),
                 currentAxes.rCut[box], ff.tolerance,
#ifdef _OPENMP
                 omp_get_max_threads());
#else
                 1);
#endif
    
    double recipPME = p.energy(molCoords.x, molCoords.y, molCoords.z,
                               atomForceRec.x, atomForceRec.y, atomForceRec.z);
    double Eextra = p.energy_extra();
    double Eself  = p.energy_self();

    // Calculate molForce based on atomForce returned from p.energy
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

    while(thisMol != end) {
      uint molIndex = *thisMol;

      uint length, start, p;
      length = mols.GetKind(molIndex).NumAtoms();
      start = mols.MolStart(molIndex);
      double X = 0.0, Y = 0.0, Z = 0.0;
      for(p = start; p < start + length; p++) {
        X += atomForceRec.x[p];
        Y += atomForceRec.y[p];
        Z += atomForceRec.z[p];
        molForceRec.Add(molIndex, X, Y, Z);
      }

      thisMol++;
    }

    return recipPME + Eextra - Eself;
  }
  return 0.0;
}

