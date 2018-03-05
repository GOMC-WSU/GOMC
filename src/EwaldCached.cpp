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

EwaldCached::EwaldCached(StaticVals & stat, System & sys) :
  forcefield(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
  currentCOM(sys.com), sysPotRef(sys.potential),
  atomForceRecRef(sys.atomForceRecRef),
  molForceRecRef(sys.molForceRecRef),
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
  multiParticleEnabled = stat.multiParticleEnabled;
}


EwaldCached::~EwaldCached()
{
  if(ewald) {
#ifdef GOMC_CUDA
    DestroyEwaldCUDAVars(forcefield.particles->getCUDAVars());
#endif
    for(int i = 0; i < mols.count; i++) {
      //when cached option is choosen
      if (cosMolRef[i] != NULL) {
        delete[] cosMolRef[i];
        delete[] sinMolRef[i];
        delete[] cosMolBoxRecip[i];
        delete[] sinMolBoxRecip[i];
      }
    }

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
      //when cached option is choosen
      if (cosMolRestore != NULL) {
        delete[] cosMolRestore;
        delete[] sinMolRestore;
      }
      //when cached option is choosen
      if (cosMolRef != NULL) {
        delete[] cosMolRef;
        delete[] sinMolRestore;
        delete[] sinMolRef;
        delete[] cosMolBoxRecip;
        delete[] sinMolBoxRecip;
      }
    }
  }
}


void EwaldCached::SetNull()
{
  //get size of image using defined Kmax
  //imageSize = imageTotal;
  //Set to NULL

  kmax = NULL;
  imageSize = NULL ;
  imageSizeRef = NULL ;
  sumRnew = NULL;
  sumInew = NULL;
  sumRref = NULL;
  sumIref = NULL;
  kx = NULL;
  ky = NULL;
  kz = NULL;
  hsqr = NULL;
  prefact = NULL;
  kxRef = NULL;
  kyRef = NULL;
  kzRef = NULL;
  hsqrRef = NULL;
  prefactRef = NULL;
  cosMolRef = NULL;
  sinMolRef = NULL;
  cosMolBoxRecip = NULL;
  sinMolBoxRecip = NULL;

  cosMolRestore = NULL;
  sinMolRestore = NULL;
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
  SetNull();
  AllocMem();
  //initialize K vectors and reciprocate terms
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    RecipInit(b, currentAxes);
    BoxReciprocalSetup(b, currentCoords);
    SetRecipRef(b);
    //calculate reciprocate force
    ForceReciprocal(atomForceRecRef, molForceRecRef, b);
    printf("Box: %d, RecipVectors: %6d, kmax: %d\n",
           b, imageSize[b], kmax[b]);
  }
}


void EwaldCached::AllocMem()
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
  cosMolRef = new double*[mols.count];
  sinMolRef = new double*[mols.count];
  cosMolBoxRecip = new double*[mols.count];
  sinMolBoxRecip = new double*[mols.count];

  for(uint b = 0; b < BOXES_WITH_U_NB; b++) {
    RecipCountInit(b, currentAxes);
  }
  //25% larger than original box size, reserved for image size change
  imageTotal = findLargeImage();
  memoryAllocation = imageTotal;

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


void EwaldCached::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  if(boxAxes.orthogonal[box])
    RecipInitOrth(box, boxAxes);
  else
    RecipInitNonOrth(box, boxAxes);
}

void EwaldCached::RecipInitOrth(uint box, BoxDimensions const& boxAxes)
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

void EwaldCached::RecipInitNonOrth(uint box, BoxDimensions const& boxAxes)
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
void EwaldCached::RecipCountInit(uint box, BoxDimensions const& boxAxes)
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


//compare number of images in different boxes and select the largest one
uint EwaldCached::findLargeImage()
{
  imageLarge = 0;
  for (int b = 0; b < BOXES_WITH_U_NB; b++) {
    if (imageLarge < imageSize[b])
      imageLarge = imageSize[b];
  }
  return imageLarge;
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

// NOTE: The calculation of W12, W13, W23 is expensive and would not be
// requied for pressure and surface tension calculation. So, they have been
// commented out. In case you need to calculate them, uncomment them.
Virial EwaldCached::VirialReciprocal(Virial& virial, uint box) const
{
  Virial tempVir = virial;
  if (box >= BOXES_WITH_U_NB)
    return tempVir;

  double wT11 = 0.0, wT12 = 0.0, wT13 = 0.0;
  double wT22 = 0.0, wT23 = 0.0, wT33 = 0.0;

  double recipIntra = 0.0;
  double constVal = 1.0 / (4.0 * alpha * alpha);
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
      double atomCharge = mols.GetKind(*thisMol).AtomCharge(p);
      chargeBox.push_back(atomCharge);
      atomIndex++;
    }
    thisMol++;
  }

  CallForceReciprocalGPU(forcefield.particles->getCUDAVars(), thisBoxCoords,
                         thisBoxCOMDiff, chargeBox, wT11, wT12,
                         wT13, wT22, wT23, wT33, imageSizeRef[box], constVal,
                         box);
#else
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(i, factor) \
  reduction(+:wT11, wT12, wT13, wT22, wT23, wT33)
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
    /*
    wT12 += factor * (-2.0 * (constVal + 1.0 / hsqrRef[box][i]) * kxRef[box][i] * kyRef[box][i]);

    wT13 += factor * (-2.0 * (constVal + 1.0 / hsqrRef[box][i]) * kxRef[box][i] * kzRef[box][i]);

    wT23 += factor * (-2.0 * (constVal + 1.0 / hsqrRef[box][i]) * kyRef[box][i] * kzRef[box][i]);
          */
  }

  //Intramolecular part
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

      // charge = particleCharge[atom];
      charge = mols.GetKind(*thisMol).AtomCharge(p);

#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(i, arg, factor) \
      reduction(+:wT11, wT12, wT13, wT22, wT23, wT33)
#endif
      for (i = 0; i < imageSizeRef[box]; i++) {
        //compute the dot product of k and r
        arg = currentAxes.DotProduct(atom, kxRef[box][i], kyRef[box][i],
                                     kzRef[box][i], currentCoords);

        factor = prefactRef[box][i] * 2.0 * (sumIref[box][i] * cos(arg) -
                                             sumRref[box][i] * sin(arg)) * charge;

        wT11 += factor * (kxRef[box][i] * diffC.x);

        wT22 += factor * (kyRef[box][i] * diffC.y);

        wT33 += factor * (kzRef[box][i] * diffC.z);

        /*
          wT12 += factor * 0.5 *(kxRef[box][i] * diffC.y + kyRef[box][i] * diffC.x);
          wT13 += factor * 0.5 *(kxRef[box][i] * diffC.z + kzRef[box][i] * diffC.x);
          wT23 += factor * 0.5 *(kyRef[box][i] * diffC.z + kzRef[box][i] * diffC.y);
        */
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


//calculate reciprocate term for displacement and rotation move
double EwaldCached::MolReciprocal(XYZArray const& molCoords,
				  const uint molIndex,
                                  const uint box)
{
  double energyRecipNew = 0.0;

  if (box < BOXES_WITH_U_NB) {
    MoleculeKind const& thisKind = mols.GetKind(molIndex);
    uint length = thisKind.NumAtoms();
    uint p, i;
    double sumRealNew, sumImaginaryNew, dotProductNew, sumRealOld,
           sumImaginaryOld;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, p, sumRealNew, sumImaginaryNew, sumRealOld, sumImaginaryOld, dotProductNew) reduction(+:energyRecipNew)
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


//calculate self term for a box
double EwaldCached::BoxSelf(BoxDimensions const& boxAxes, uint box) const
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

  self = -1.0 * self * alpha * num::qqFact / sqrt(M_PI);

  return self;
}


//calculate correction term for a molecule
double EwaldCached::MolCorrection(uint molIndex, uint box) const
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
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords,
                         start + i, start + j, box);
      dist = sqrt(distSq);
      correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                     erf(alpha * dist) / dist);
    }
  }

  return correction;
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


double EwaldCached::SwapSelf(const cbmc::TrialMol& trialMol) const
{
  if (trialMol.GetBox() >= BOXES_WITH_U_NB)
    return 0.0;

  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  double en_self = 0.0;

  for (uint i = 0; i < atomSize; i++) {
    en_self -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(i));
  }
  return (en_self * alpha * num::qqFact / sqrt(M_PI));
}

//calculate correction term for linear molecule CBMC algorithm
double EwaldCached::SwapCorrection(const cbmc::TrialMol& trialMol) const
{
  if (trialMol.GetBox() >= BOXES_WITH_U_NB)
    return 0.0;

  double dist, distSq;
  double correction = 0.0;
  XYZ virComponents;
  const MoleculeKind& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();

  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, trialMol.GetBox());

      dist = sqrt(distSq);
      correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                     erf(alpha * dist) / dist);
    }
  }
  return num::qqFact * correction;
}


//back up reciptocate value to Ref (will be called during initialization)
void EwaldCached::SetRecipRef(uint box)
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


//update reciprocate values
void EwaldCached::UpdateRecip(uint box)
{
  double *tempR, *tempI;
  tempR = sumRref[box];
  tempI = sumIref[box];
  sumRref[box] = sumRnew[box];
  sumIref[box] = sumInew[box];
  sumRnew[box] = tempR;
  sumInew[box] = tempI;
#ifdef GOMC_CUDA
  UpdateRecipCUDA(forcefield.particles->getCUDAVars(), box);
#endif
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


//calculate reciprocate force term for a box with molCoords
void EwaldCached::BoxForceReciprocal(XYZArray const& molCoords,
                                     XYZArray& atomForceRec,
                                     XYZArray& molForceRec,
                                     uint box)
{
  if(multiParticleEnabled && (box < BOXES_WITH_U_NB)) {
    // molecule iterator
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

    while(thisMol != end) {
      uint molIndex = *thisMol;
      uint length, start, p, i;
      double dot, factor;
      molForceRec.Set(molIndex, 0.0, 0.0, 0.0);
      length = mols.GetKind(molIndex).NumAtoms();
      start = mols.MolStart(molIndex);

      for(p = start; p < start + length; p++) {
        double X = 0.0, Y = 0.0, Z = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(i, dot, factor) \
        reduction(+:X, Y, Z)
#endif
        for(i = 0; i < imageSize[box]; i++) {
          dot = currentAxes.DotProduct(p, kx[box][i], ky[box][i], kz[box][i],
                                       molCoords);
  
          factor = particleCharge[p] * prefact[box][i] *
                  (sin(dot) * sumRnew[box][i] - cos(dot) * sumInew[box][i]);

          X += factor * kx[box][i];
          Y += factor * ky[box][i];
          Z += factor * kz[box][i];
        }
        atomForceRec.Set(p, X, Y, Z);
        molForceRec.Add(molIndex, X, Y, Z);
      }
      thisMol++;
    }
  }
}

//calculate reciprocate force term for a box with Reference value
void EwaldCached::ForceReciprocal(XYZArray& atomForceRec, XYZArray& molForceRec,
				                          uint box)
{
  if(multiParticleEnabled && (box < BOXES_WITH_U_NB)) {
    // molecule iterator
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

    while(thisMol != end) {
      uint molIndex = *thisMol;
      uint length, start, p, i;
      double dot, factor;
      length = mols.GetKind(molIndex).NumAtoms();
      start = mols.MolStart(molIndex);
      molForceRec.Set(molIndex, 0.0, 0.0, 0.0);

      for(p = start; p < start + length; p++) {
        double X = 0.0, Y = 0.0, Z = 0.0;
#ifdef _OPENMP
        #pragma omp parallel for default(shared) \
        private(i, dot, factor) reduction(+:X, Y, Z)
#endif
        for(i = 0; i < imageSizeRef[box]; i++) {
          dot = currentAxes.DotProduct(p, kxRef[box][i], kyRef[box][i],
                                        kzRef[box][i], currentCoords);
    
          factor = particleCharge[p] * prefactRef[box][i] *
                    (sin(dot) * sumRref[box][i] - cos(dot) * sumIref[box][i]);
    
          X += factor * kxRef[box][i];
          Y += factor * kyRef[box][i];
          Z += factor * kzRef[box][i];
        }
        atomForceRec.Set(p, X, Y, Z);
        molForceRec.Add(molIndex, X, Y, Z);
      }
      thisMol++;
    }
  }
}


