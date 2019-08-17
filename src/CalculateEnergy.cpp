/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "CalculateEnergy.h"        //header for this
#include "EwaldCached.h"            //for ewald calculation
#include "Ewald.h"                  //for ewald calculation
#include "NoEwald.h"                //for ewald calculation
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
#include "BoxDimensionsNonOrth.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#include <algorithm>
#ifdef GOMC_CUDA
#include "CalculateEnergyCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

#define HALF_DOUBLE_MAX 8.988466e+307

using namespace geom;

CalculateEnergy::CalculateEnergy(StaticVals & stat, System & sys) :
  forcefield(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
  currentCOM(sys.com),
  atomForceRef(sys.atomForceRef),
  molForceRef(sys.molForceRef),
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
  , cellList(sys.cellList)
{
}


void CalculateEnergy::Init(System & sys)
{
  uint maxAtomInMol = 0;
  calcEwald = sys.GetEwald();
  electrostatic = forcefield.electrostatic;
  ewald = forcefield.ewald;
  multiParticleEnabled = sys.statV.multiParticleEnabled;
  for(uint m = 0; m < mols.count; ++m) {
    const MoleculeKind& molKind = mols.GetKind(m);
    if(molKind.NumAtoms() > maxAtomInMol)
      maxAtomInMol = molKind.NumAtoms();
    for(uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleKind.push_back(molKind.AtomKind(a));
      particleMol.push_back(m);
      particleCharge.push_back(molKind.AtomCharge(a));
    }
  }
  // initialize the energy table with pre-calculated energies
  if(sys.statV.forcefield.energyTable) {
    energyTableEnabled = true;
    forcefield.particles->InitializeTables();
  } else {
    energyTableEnabled = false;
  }
#ifdef GOMC_CUDA
  InitCoordinatesCUDA(forcefield.particles->getCUDAVars(),
                      currentCoords.Count(), maxAtomInMol, currentCOM.Count());
#endif
}

SystemPotential CalculateEnergy::SystemTotal()
{
  SystemPotential pot =
    SystemInter(SystemPotential(), currentCoords, currentCOM, atomForceRef,
                molForceRef, currentAxes);

  //system intra
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    int i;
    double bondEnergy[2] = {0};
    double bondEn = 0.0, nonbondEn = 0.0, self = 0.0, correction = 0.0;
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(b);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(b);
    std::vector<uint> molID;

    while (thisMol != end) {
      molID.push_back(*thisMol);
      ++thisMol;
    }

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, bondEnergy) reduction(+:bondEn, nonbondEn, correction)
#endif
    for (i = 0; i < molID.size(); i++) {
      //calculate nonbonded energy
      MoleculeIntra(molID[i], b, bondEnergy);
      bondEn += bondEnergy[0];
      nonbondEn += bondEnergy[1];
      //calculate correction term of electrostatic interaction
      correction += calcEwald->MolCorrection(molID[i], b);
    }

    pot.boxEnergy[b].intraBond = bondEn;
    pot.boxEnergy[b].intraNonbond = nonbondEn;
    //calculate self term of electrostatic interaction
    pot.boxEnergy[b].self = calcEwald->BoxSelf(currentAxes, b);
    pot.boxEnergy[b].correction = -1 * correction * num::qqFact;

    //Calculate Virial
    pot.boxVirial[b] = VirialCalc(b);
  }

  pot.Total();

  if(pot.totalEnergy.total > 1.0e12) {
    std::cout << "\nWarning: Large energy detected due to the overlap in "
              "initial configuration.\n"
              "         The total energy will be recalculated at EqStep to "  
              "ensure the accuracy \n"
              "         of the computed running energies.\n";
  }

  return pot;
}


SystemPotential CalculateEnergy::SystemInter(SystemPotential potential,
  XYZArray const& coords,
  XYZArray const& com,
  XYZArray& atomForce,
  XYZArray& molForce,
  BoxDimensions const& boxAxes)
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    //calculate LJ interaction and real term of electrostatic interaction
    potential = BoxInter(potential, coords, atomForce, molForce,
                         boxAxes, b);
    //calculate reciprocate term of electrostatic interaction
    potential.boxEnergy[b].recip = calcEwald->BoxReciprocal(b);
  }

  potential.Total();

  return potential;
}

SystemPotential CalculateEnergy::BoxInter(SystemPotential potential,
                                          XYZArray const& coords,
                                          XYZArray& atomForce,
                                          XYZArray& molForce,
                                          BoxDimensions const& boxAxes,
                                          const uint box)
{
  //Handles reservoir box case, returning zeroed structure if
  //interactions are off.
  if (box >= BOXES_WITH_U_NB)
    return potential;

  double tempREn = 0.0, tempLJEn = 0.0;
  double distSq, qi_qj_fact;
  int i;
  XYZ virComponents, force, forceReal, forceLJ;
  std::vector<uint> pair1, pair2;
  CellList::Pairs pair = cellList.EnumeratePairs(box);
  double rREn, rLJEn, rCoulomb, rForce;

  // make a pointer to atom force and mol force for openmp
  double *aForcex = atomForce.x;
  double *aForcey = atomForce.y;
  double *aForcez = atomForce.z;
  double *mForcex = molForce.x;
  double *mForcey = molForce.y;
  double *mForcez = molForce.z;
  int atomCount = atomForce.Count();
  int molCount = molForce.Count();

  // Reset Force Arrays
  ResetForce(atomForce, molForce, box);

  //store atom pair index
  while (!pair.Done()) {
    if(!SameMolecule(pair.First(), pair.Second())) {
      pair1.push_back(pair.First());
      pair2.push_back(pair.Second());
    }
    pair.Next();
  }

#ifdef GOMC_CUDA
  uint pairSize = pair1.size();
  uint currentIndex = 0;
  double REn = 0.0, LJEn = 0.0;
  //update unitcell in GPU
  UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                      boxAxes.cellBasis[box].x, boxAxes.cellBasis[box].y,
                      boxAxes.cellBasis[box].z);

  if(!boxAxes.orthogonal[box]) {
    BoxDimensionsNonOrth newAxes = *((BoxDimensionsNonOrth*)(&boxAxes));
    UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                           newAxes.cellBasis_Inv[box].x,
                           newAxes.cellBasis_Inv[box].y,
                           newAxes.cellBasis_Inv[box].z);
  }

  while(currentIndex < pairSize) {
    uint max = currentIndex + MAX_PAIR_SIZE;
    max = (max < pairSize ? max : pairSize);

    std::vector<uint>::const_iterator first1 = pair1.begin() + currentIndex;
    std::vector<uint>::const_iterator last1 = pair1.begin() + max;
    std::vector<uint>::const_iterator first2 = pair2.begin() + currentIndex;
    std::vector<uint>::const_iterator last2 = pair2.begin() + max;
    std::vector<uint> subPair1(first1, last1);
    std::vector<uint> subPair2(first2, last2);

    // Reset forces on GPU for the first iteration
    bool reset_force = currentIndex == 0;

    // Copy back the result if it is the last iteration
    bool copy_back = max == pairSize;

    CallBoxInterGPU(forcefield.particles->getCUDAVars(), subPair1, subPair2,
                    coords, boxAxes, electrostatic, particleCharge,
                    particleKind, particleMol, REn, LJEn, multiParticleEnabled,
                    aForcex, aForcey, aForcez, mForcex, mForcey, mForcez,
                    atomCount, molCount, reset_force, copy_back, box);
    tempREn += REn;
    tempLJEn += LJEn;
    currentIndex += MAX_PAIR_SIZE;
  }

#else
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, virComponents, forceReal, forceLJ) \
reduction(+:tempREn, tempLJEn, \
aForcex[:atomCount], aForcey[:atomCount], aForcez[:atomCount], \
mForcex[:molCount], mForcey[:molCount], mForcez[:molCount])
#endif
  for (i = 0; i < pair1.size(); i++) {
    if(boxAxes.InRcut(distSq, virComponents, coords, pair1[i], pair2[i], box)) {
      if (electrostatic) {
        qi_qj_fact = particleCharge[pair1[i]] * particleCharge[pair2[i]] * num::qqFact;
        tempREn += forcefield.particles->CalcCoulomb(distSq, qi_qj_fact, box);
      }
      tempLJEn += forcefield.particles->CalcEn(distSq, particleKind[pair1[i]], particleKind[pair2[i]]);

      // In case of multiparticle we also need to calculate force
      if(multiParticleEnabled) {
        if(electrostatic) {
          forceReal = virComponents *
            forcefield.particles->CalcCoulombVir(distSq, qi_qj_fact, box);
        }
        forceLJ = virComponents *
          forcefield.particles->CalcVir(distSq, particleKind[pair1[i]],
                                        particleKind[pair2[i]]);
        aForcex[pair1[i]] += forceLJ.x + forceReal.x;
        aForcey[pair1[i]] += forceLJ.y + forceReal.y;
        aForcez[pair1[i]] += forceLJ.z + forceReal.z;
        aForcex[pair2[i]] += -(forceLJ.x + forceReal.x);
        aForcey[pair2[i]] += -(forceLJ.y + forceReal.y);
        aForcez[pair2[i]] += -(forceLJ.z + forceReal.z);
        mForcex[particleMol[pair1[i]]] += (forceLJ.x + forceReal.x);
        mForcey[particleMol[pair1[i]]] += (forceLJ.y + forceReal.y);
        mForcez[particleMol[pair1[i]]] += (forceLJ.z + forceReal.z);
        mForcex[particleMol[pair2[i]]] += -(forceLJ.x + forceReal.x);
        mForcey[particleMol[pair2[i]]] += -(forceLJ.y + forceReal.y);
        mForcez[particleMol[pair2[i]]] += -(forceLJ.z + forceReal.z);
      }
    }
  }
#endif

  // setting energy and virial of LJ interaction
  potential.boxEnergy[box].inter = tempLJEn;
  // setting energy and virial of coulomb interaction
  potential.boxEnergy[box].real = tempREn;

  // set correction energy and virial
  if (forcefield.useLRC) {
    EnergyCorrection(potential, boxAxes, box);
  }

  potential.Total();

  return potential;
}

// NOTE: The calculation of W12, W13, W23 is expensive and would not be
// requied for pressure and surface tension calculation. So, they have been
// commented out. In case you need to calculate them, uncomment them.
Virial CalculateEnergy::VirialCalc(const uint box)
{
  //store virial and energy of reference and modify the virial
  Virial tempVir;

  //tensors for VDW and real part of electrostatic
  double vT11 = 0.0, vT12 = 0.0, vT13 = 0.0;
  double vT22 = 0.0, vT23 = 0.0, vT33 = 0.0;
  double rT11 = 0.0, rT12 = 0.0, rT13 = 0.0;
  double rT22 = 0.0, rT23 = 0.0, rT33 = 0.0;

  double distSq, pVF, pRF, qi_qj;
  int i;
  XYZ virC, comC;
  std::vector<uint> pair1, pair2;
  CellList::Pairs pair = cellList.EnumeratePairs(box);
  //store atom pair index
  while (!pair.Done()) {
    if(!SameMolecule(pair.First(), pair.Second())) {
      pair1.push_back(pair.First());
      pair2.push_back(pair.Second());
    }
    pair.Next();
  }

#ifdef GOMC_CUDA
  uint pairSize = pair1.size();
  //update unitcell in GPU
  UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                      currentAxes.cellBasis[box].x,
                      currentAxes.cellBasis[box].y,
                      currentAxes.cellBasis[box].z);

  if(!currentAxes.orthogonal[box]) {
    BoxDimensionsNonOrth newAxes = *((BoxDimensionsNonOrth*)(&currentAxes));
    UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                           newAxes.cellBasis_Inv[box].x, newAxes.cellBasis_Inv[box].y,
                           newAxes.cellBasis_Inv[box].z);
  }

  uint currentIndex = 0;
  double vT11t = 0.0, vT12t = 0.0, vT13t = 0.0;
  double vT22t = 0.0, vT23t = 0.0, vT33t = 0.0;
  double rT11t = 0.0, rT12t = 0.0, rT13t = 0.0;
  double rT22t = 0.0, rT23t = 0.0, rT33t = 0.0;
  while(currentIndex < pairSize) {
    uint max = currentIndex + MAX_PAIR_SIZE;
    max = (max < pairSize ? max : pairSize);

    std::vector<uint>::const_iterator first1 = pair1.begin() + currentIndex;
    std::vector<uint>::const_iterator last1 = pair1.begin() + max;
    std::vector<uint>::const_iterator first2 = pair2.begin() + currentIndex;
    std::vector<uint>::const_iterator last2 = pair2.begin() + max;
    std::vector<uint> subPair1(first1, last1);
    std::vector<uint> subPair2(first2, last2);
    CallBoxInterForceGPU(forcefield.particles->getCUDAVars(), subPair1,
                         subPair2, currentCoords, currentCOM, currentAxes,
                         electrostatic, particleCharge, particleKind,
                         particleMol, rT11t, rT12t, rT13t, rT22t, rT23t, rT33t,
                         vT11t, vT12t, vT13t, vT22t, vT23t, vT33t, box);
    rT11 += rT11t;
    rT12 += rT12t;
    rT13 += rT13t;
    rT22 += rT22t;
    rT23 += rT23t;
    rT33 += rT33t;
    vT11 += vT11t;
    vT12 += vT12t;
    vT13 += vT13t;
    vT22 += vT22t;
    vT23 += vT23t;
    vT33 += vT33t;
    currentIndex += MAX_PAIR_SIZE;
  }
#else
#ifdef _OPENMP
  #pragma omp parallel for default(shared) private(i, distSq, pVF, pRF, qi_qj, virC, comC) reduction(+:vT11, vT12, vT13, vT22, vT23, vT33, rT11, rT12, rT13, rT22, rT23, rT33)
#endif
  for (i = 0; i < pair1.size(); i++) {
    if (currentAxes.InRcut(distSq, virC, currentCoords, pair1[i],
                           pair2[i], box)) {
      pVF = 0.0;
      pRF = 0.0;

      //calculate the distance between com of two molecules
      comC = currentCOM.Difference(particleMol[pair1[i]],
                                   particleMol[pair2[i]]);
      //calculate the minimum image between com of two molecules
      comC = currentAxes.MinImage(comC, box);

      if (electrostatic) {
        qi_qj = particleCharge[pair1[i]] * particleCharge[pair2[i]];

        pRF = forcefield.particles->CalcCoulombVir(distSq, qi_qj, box);
        //calculate the top diagonal of pressure tensor
        rT11 += pRF * (virC.x * comC.x);
        //rT12 += pRF * (0.5 * (virC.x * comC.y + virC.y * comC.x));
        //rT13 += pRF * (0.5 * (virC.x * comC.z + virC.z * comC.x));

        rT22 += pRF * (virC.y * comC.y);
        //rT23 += pRF * (0.5 * (virC.y * comC.z + virC.z * comC.y));

        rT33 += pRF * (virC.z * comC.z);
      }

      pVF = forcefield.particles->CalcVir(distSq, particleKind[pair1[i]],
                                          particleKind[pair2[i]]);
      //calculate the top diagonal of pressure tensor
      vT11 += pVF * (virC.x * comC.x);
      //vT12 += pVF * (0.5 * (virC.x * comC.y + virC.y * comC.x));
      //vT13 += pVF * (0.5 * (virC.x * comC.z + virC.z * comC.x));

      vT22 += pVF * (virC.y * comC.y);
      //vT23 += pVF * (0.5 * (virC.y * comC.z + virC.z * comC.y));

      vT33 += pVF * (virC.z * comC.z);
    }
  }
#endif

  // set the all tensor values
  tempVir.interTens[0][0] = vT11;
  tempVir.interTens[0][1] = vT12;
  tempVir.interTens[0][2] = vT13;

  tempVir.interTens[1][0] = vT12;
  tempVir.interTens[1][1] = vT22;
  tempVir.interTens[1][2] = vT23;

  tempVir.interTens[2][0] = vT13;
  tempVir.interTens[2][1] = vT23;
  tempVir.interTens[2][2] = vT33;

  if (electrostatic) {
    // real part of electrostatic
    tempVir.realTens[0][0] = rT11 * num::qqFact;
    tempVir.realTens[0][1] = rT12 * num::qqFact;
    tempVir.realTens[0][2] = rT13 * num::qqFact;

    tempVir.realTens[1][0] = rT12 * num::qqFact;
    tempVir.realTens[1][1] = rT22 * num::qqFact;
    tempVir.realTens[1][2] = rT23 * num::qqFact;

    tempVir.realTens[2][0] = rT13 * num::qqFact;
    tempVir.realTens[2][1] = rT23 * num::qqFact;
    tempVir.realTens[2][2] = rT33 * num::qqFact;
  }

  // setting virial of LJ
  tempVir.inter = vT11 + vT22 + vT33;
  // setting virial of coulomb
  tempVir.real = (rT11 + rT22 + rT33) * num::qqFact;

  if (forcefield.useLRC) {
    VirialCorrection(tempVir, currentAxes, box);
  }

  //calculate reciprocate term of force
  tempVir = calcEwald->VirialReciprocal(tempVir, box);

  tempVir.Total();

  return tempVir;
}

bool CalculateEnergy::MoleculeInter(Intermolecular &inter_LJ,
                                    Intermolecular &inter_coulomb,
                                    XYZArray const& molCoords,
                                    const uint molIndex,
                                    const uint box) const
{
  double tempREn = 0.0, tempLJEn = 0.0;
  bool overlap = false;

  if (box < BOXES_WITH_U_NB) {
    uint length = mols.GetKind(molIndex).NumAtoms();
    uint start = mols.MolStart(molIndex);

    for (uint p = 0; p < length; ++p) {
      uint atom = start + p;
      CellList::Neighbors n = cellList.EnumerateLocal(currentCoords[atom],
                              box);
      n = cellList.EnumerateLocal(currentCoords[atom], box);

      double qi_qj_fact, distSq;
      int i;
      XYZ virComponents, forceLJ, forceReal;
      std::vector<uint> nIndex;

      //store atom index in neighboring cell
      while (!n.Done()) {
        nIndex.push_back(*n);
        n.Next();
      }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, virComponents, forceLJ, forceReal) reduction(+:tempREn, tempLJEn)
#endif
      for(i = 0; i < nIndex.size(); i++) {
        distSq = 0.0;
        //Subtract old energy
        if (currentAxes.InRcut(distSq, virComponents,
                               currentCoords, atom, nIndex[i], box)) {

          if (electrostatic) {
            qi_qj_fact = particleCharge[atom] * particleCharge[nIndex[i]] *
                         num::qqFact;

            tempREn -= forcefield.particles->CalcCoulomb(distSq,qi_qj_fact,box);
          }

          tempLJEn -= forcefield.particles->CalcEn(distSq, particleKind[atom],
                      particleKind[nIndex[i]]);
        }
      }

      //add new energy
      n = cellList.EnumerateLocal(molCoords[p], box);
      //store atom index in neighboring cell
      nIndex.clear();
      while (!n.Done()) {
        nIndex.push_back(*n);
        n.Next();
      }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, virComponents, forceReal, forceLJ) reduction(+:tempREn, tempLJEn)
#endif
      for(i = 0; i < nIndex.size(); i++) {
        distSq = 0.0;
        if (currentAxes.InRcut(distSq, virComponents,
                               molCoords, p, currentCoords, nIndex[i], box)) {
          if(distSq < forcefield.rCutLowSq) {
            overlap |= true;
          }

          if (electrostatic) {
            qi_qj_fact = particleCharge[atom] *
                         particleCharge[nIndex[i]] * num::qqFact;

            tempREn += forcefield.particles->CalcCoulomb(distSq,
                       qi_qj_fact, box);
          }

          tempLJEn += forcefield.particles->CalcEn(distSq,
                      particleKind[atom],
                      particleKind[nIndex[i]]);
        }
      }
    }
  }

  inter_LJ.energy = tempLJEn;
  inter_coulomb.energy = tempREn;
  return overlap;
}

// Calculate 1-N nonbonded intra energy
void CalculateEnergy::ParticleNonbonded(double* inter,
                                        cbmc::TrialMol const& trialMol,
                                        XYZArray const& trialPos,
                                        const uint partIndex,
                                        const uint box,
                                        const uint trials) const
{
  if (box >= BOXES_WITH_U_B)
    return;

  const MoleculeKind& kind = trialMol.GetKind();
  //loop over all partners of the trial particle
  const uint* partner = kind.sortedNB.Begin(partIndex);
  const uint* end = kind.sortedNB.End(partIndex);
  while (partner != end) {
    if (trialMol.AtomExists(*partner)) {
      for (uint t = 0; t < trials; ++t) {
        double distSq;

        if (currentAxes.InRcut(distSq, trialPos, t, trialMol.GetCoords(),
                               *partner, box)) {                     
          inter[t] += forcefield.particles->CalcEn(distSq,
                      kind.AtomKind(partIndex),
                      kind.AtomKind(*partner));
          if (electrostatic) {
            double qi_qj_Fact = kind.AtomCharge(partIndex) *
                                kind.AtomCharge(*partner) * num::qqFact;
            forcefield.particles->CalcCoulombAdd_1_4(inter[t], distSq,
                qi_qj_Fact, true);
          }
        }
      }
    }
    ++partner;
  }
}

void CalculateEnergy::ParticleInter(double* en, double *real,
                                    XYZArray const& trialPos,
                                    bool* overlap,
                                    const uint partIndex,
                                    const uint molIndex,
                                    const uint box,
                                    const uint trials) const
{
  if(box >= BOXES_WITH_U_NB)
    return;
  double distSq, qi_qj_Fact, tempLJ, tempReal;
  int i;
  MoleculeKind const& thisKind = mols.GetKind(molIndex);
  uint kindI = thisKind.AtomKind(partIndex);
  double kindICharge = thisKind.AtomCharge(partIndex);
  std::vector<uint> nIndex;

  for(uint t = 0; t < trials; ++t) {
    nIndex.clear();
    tempReal = 0.0;
    tempLJ = 0.0;
    CellList::Neighbors n = cellList.EnumerateLocal(trialPos[t], box);
    while (!n.Done()) {
      nIndex.push_back(*n);
      n.Next();
    }

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(i, distSq, qi_qj_Fact) reduction(+:tempLJ, tempReal)
#endif
    for(i = 0; i < nIndex.size(); i++) {
      distSq = 0.0;

      if(currentAxes.InRcut(distSq, trialPos, t, currentCoords, nIndex[i],box)) {
        if(distSq < forcefield.rCutLowSq) {
          overlap[t] |= true;
        } 
        double tt = forcefield.particles->CalcEn(distSq, kindI,
                                               particleKind[nIndex[i]]);
        if(energyTableEnabled && isnan(tt)) {
          // In energy table we lose some accuracy, so some energies for very
          // small distSq can be nan. In those cases we set them to a very large
          // number so the program continues as normal
          tt = HALF_DOUBLE_MAX;
        }
        tempLJ += tt;
        if(electrostatic) {
          qi_qj_Fact = particleCharge[nIndex[i]] * kindICharge * num::qqFact;
          tempReal += forcefield.particles->CalcCoulomb(distSq, qi_qj_Fact, box);
        }
      }
    }
    en[t] += tempLJ;
    real[t] += tempReal;
  }
}


//Calculates the change in the TC from adding numChange atoms of a kind
Intermolecular CalculateEnergy::MoleculeTailChange(const uint box,
                                                   const uint kind,
                                                   const bool add) const
{
  Intermolecular delta;

  if (box < BOXES_WITH_U_NB) {

    double sign = (add ? 1.0 : -1.0);
    uint mkIdxII = kind * mols.GetKindsCount() + kind;
    for (uint j = 0; j < mols.GetKindsCount(); ++j) {
      uint mkIdxIJ = j * mols.GetKindsCount() + kind;
      double rhoDeltaIJ_2 = sign * 2.0 *
                            (double)(molLookup.NumKindInBox(j, box)) * currentAxes.volInv[box];
      delta.energy += mols.pairEnCorrections[mkIdxIJ] * rhoDeltaIJ_2;
    }

    //We already calculated part of the change for this type in the loop
    delta.energy += mols.pairEnCorrections[mkIdxII] *
                    currentAxes.volInv[box];
  }
  return delta;
}


//Calculates intramolecular energy of a full molecule
void CalculateEnergy::MoleculeIntra(const uint molIndex,
                                    const uint box, double *bondEn) const
{
  bondEn[0] = 0.0, bondEn[1] = 0.0;

  MoleculeKind& molKind = mols.kinds[mols.kIndex[molIndex]];
  // *2 because we'll be storing inverse bond vectors
  XYZArray bondVec(molKind.bondList.count * 2);

  BondVectors(bondVec, molKind, molIndex, box);
  MolBond(bondEn[0], molKind, bondVec, molIndex, box);
  MolAngle(bondEn[0], molKind, bondVec, box);
  MolDihedral(bondEn[0], molKind, bondVec, box);
  MolNonbond(bondEn[1], molKind, molIndex, box);
  MolNonbond_1_4(bondEn[1], molKind, molIndex, box);
  MolNonbond_1_3(bondEn[1], molKind, molIndex, box);
}

//used in molecule exchange for calculating bonded and intraNonbonded energy
Energy CalculateEnergy::MoleculeIntra(cbmc::TrialMol const &mol,
                                      const uint molIndex) const
{
  double bondEn = 0.0, intraNonbondEn = 0.0;
  // *2 because we'll be storing inverse bond vectors
  const MoleculeKind& molKind = mol.GetKind();
  uint count = molKind.bondList.count;
  XYZArray bondVec(count * 2);
  std::vector<bool> bondExist(count * 2, false);

  BondVectors(bondVec, mol, bondExist, molKind);  
  MolBond(bondEn, mol, bondVec, bondExist, molKind);
  MolAngle(bondEn, mol, bondVec, bondExist, molKind);
  MolDihedral(bondEn, mol, bondVec, bondExist, molKind);
  MolNonbond(intraNonbondEn, mol, molKind);
  MolNonbond_1_4(intraNonbondEn, mol, molKind);
  MolNonbond_1_3(intraNonbondEn, mol, molKind);
  return Energy(bondEn, intraNonbondEn, 0.0, 0.0, 0.0, 0.0, 0.0);
}

void CalculateEnergy::BondVectors(XYZArray & vecs,
                                  MoleculeKind const& molKind,
                                  const uint molIndex,
                                  const uint box) const
{
  for (uint i = 0; i < molKind.bondList.count; ++i) {
    uint p1 = mols.start[molIndex] + molKind.bondList.part1[i];
    uint p2 = mols.start[molIndex] + molKind.bondList.part2[i];
    XYZ dist = currentCoords.Difference(p2, p1);
    dist = currentAxes.MinImage(dist, box);

    //store inverse vectors at i+count
    vecs.Set(i, dist);
    vecs.Set(i + molKind.bondList.count, -dist.x, -dist.y, -dist.z);
  }
}

void CalculateEnergy::BondVectors(XYZArray & vecs,
                                  cbmc::TrialMol const &mol,
                                  std::vector<bool> & bondExist,
                                  MoleculeKind const& molKind) const
{
  uint box = mol.GetBox();
  uint count = molKind.bondList.count;
  for (uint i = 0; i < count; ++i) {
    uint p1 = molKind.bondList.part1[i];
    uint p2 = molKind.bondList.part2[i];
    if(mol.AtomExists(p1) && mol.AtomExists(p2)) {
      bondExist[i] = true;
      bondExist[i + count] = true;
      XYZ dist = mol.GetCoords().Difference(p2, p1);
      dist = currentAxes.MinImage(dist, box);
      //store inverse vectors at i+count
      vecs.Set(i, dist);
      vecs.Set(i + count, -dist.x, -dist.y, -dist.z);
    }
  }
}


void CalculateEnergy::MolBond(double & energy,
                              MoleculeKind const& molKind,
                              XYZArray const& vecs,
                              const uint molIndex,
                              const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;

  for (uint b = 0; b < molKind.bondList.count; ++b) {
    double molLength = vecs.Get(b).Length();
    double eqLength = forcefield.bonds.Length(molKind.bondList.kinds[b]);
    energy += forcefield.bonds.Calc(molKind.bondList.kinds[b], molLength);
    if(abs(molLength - eqLength) > 0.02) {
      uint p1 = molKind.bondList.part1[b];
      uint p2 = molKind.bondList.part2[b];
      printf("Warning: Box%d, %6d %4s,", box, molIndex, molKind.name.c_str());
      printf("%3s-%-3s bond: Par-file ", molKind.atomNames[p1].c_str(),
            molKind.atomNames[p2].c_str());
      printf("%2.3f A, PDB file %2.3f A!\n", eqLength, molLength);
    }
  }
}

void CalculateEnergy::MolBond(double & energy,
                              cbmc::TrialMol const &mol,
                              XYZArray const& vecs,
                              std::vector<bool> const & bondExist,
                              MoleculeKind const& molKind)const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;

  uint count = molKind.bondList.count;
  for (uint b = 0; b < count; ++b) {
    if(bondExist[b]) {
      energy += forcefield.bonds.Calc(molKind.bondList.kinds[b],
                                      vecs.Get(b).Length());
    }
  }
}

void CalculateEnergy::MolAngle(double & energy,
                               MoleculeKind const& molKind,
                               XYZArray const& vecs,
                               const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;
  for (uint a = 0; a < molKind.angles.Count(); ++a) {
    //Note: need to reverse the second bond to get angle properly.
    double theta = Theta(vecs.Get(molKind.angles.GetBond(a, 0)),
                         -vecs.Get(molKind.angles.GetBond(a, 1)));
    energy += forcefield.angles->Calc(molKind.angles.GetKind(a), theta);
  }
}

void CalculateEnergy::MolAngle(double & energy,
                               cbmc::TrialMol const &mol,
                               XYZArray const& vecs,
                               std::vector<bool> const & bondExist,
                               MoleculeKind const& molKind) const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;

  uint count = molKind.angles.Count();
  for (uint a = 0; a < count; ++a) {
    if(bondExist[molKind.angles.GetBond(a, 0)] &&
      bondExist[molKind.angles.GetBond(a, 1)]) {
      //Note: need to reverse the second bond to get angle properly.
      double theta = Theta(vecs.Get(molKind.angles.GetBond(a, 0)),
                          -vecs.Get(molKind.angles.GetBond(a, 1)));
      energy += forcefield.angles->Calc(molKind.angles.GetKind(a), theta);
    }
  }
}

void CalculateEnergy::MolDihedral(double & energy,
                                  MoleculeKind const& molKind,
                                  XYZArray const& vecs,
                                  const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;
  for (uint d = 0; d < molKind.dihedrals.Count(); ++d) {
    double phi = Phi(vecs.Get(molKind.dihedrals.GetBond(d, 0)),
                     vecs.Get(molKind.dihedrals.GetBond(d, 1)),
                     vecs.Get(molKind.dihedrals.GetBond(d, 2)));
    energy += forcefield.dihedrals.Calc(molKind.dihedrals.GetKind(d), phi);
  }
}

void CalculateEnergy::MolDihedral(double & energy,
                                  cbmc::TrialMol const &mol,
                                  XYZArray const& vecs,
                                  std::vector<bool> const & bondExist,
                                  MoleculeKind const& molKind) const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;
  
  uint count =  molKind.dihedrals.Count();
  for (uint d = 0; d < count; ++d) {
    if(bondExist[molKind.dihedrals.GetBond(d, 0)] &&
      bondExist[molKind.dihedrals.GetBond(d, 1)] &&
      bondExist[molKind.dihedrals.GetBond(d, 2)]) {
      double phi = Phi(vecs.Get(molKind.dihedrals.GetBond(d, 0)),
                      vecs.Get(molKind.dihedrals.GetBond(d, 1)),
                      vecs.Get(molKind.dihedrals.GetBond(d, 2)));
      energy += forcefield.dihedrals.Calc(molKind.dihedrals.GetKind(d), phi);
    }
  }
}

// Calculate 1-N nonbonded intra energy
void CalculateEnergy::MolNonbond(double & energy,
                                 MoleculeKind const& molKind,
                                 const uint molIndex,
                                 const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;

  for (uint i = 0; i < molKind.nonBonded.count; ++i) {
    uint p1 = mols.start[molIndex] + molKind.nonBonded.part1[i];
    uint p2 = mols.start[molIndex] + molKind.nonBonded.part2[i];
    currentAxes.InRcut(distSq, currentCoords, p1, p2, box);
    if (forcefield.rCutSq > distSq) {
      energy += forcefield.particles->CalcEn(distSq, molKind.AtomKind
                                             (molKind.nonBonded.part1[i]),
                                             molKind.AtomKind
                                             (molKind.nonBonded.part2[i]));
      if (electrostatic) {
        qi_qj_Fact = num::qqFact *
                     molKind.AtomCharge(molKind.nonBonded.part1[i]) *
                     molKind.AtomCharge(molKind.nonBonded.part2[i]);

        forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
            qi_qj_Fact, true);
      }
    }
  }

}

// Calculate 1-N nonbonded intra energy uisng pos
void CalculateEnergy::MolNonbond(double & energy, cbmc::TrialMol const &mol,
                                 MoleculeKind const& molKind) const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;
  uint count = molKind.nonBonded.count;

  for (uint i = 0; i < count; ++i) {
    uint p1 = molKind.nonBonded.part1[i];
    uint p2 = molKind.nonBonded.part2[i];
    if(mol.AtomExists(p1) && mol.AtomExists(p2)) {
      currentAxes.InRcut(distSq, mol.GetCoords(), p1, p2, mol.GetBox());
      if (forcefield.rCutSq > distSq) {
        energy += forcefield.particles->CalcEn(distSq, molKind.AtomKind(p1),
                                               molKind.AtomKind(p2));
        if (electrostatic) {
          qi_qj_Fact = num::qqFact * molKind.AtomCharge(1) *
                      molKind.AtomCharge(p2);

          forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
              qi_qj_Fact, true);
        }
      }
    }
  }

}

// Calculate 1-4 nonbonded intra energy
void CalculateEnergy::MolNonbond_1_4(double & energy,
                                     MoleculeKind const& molKind,
                                     const uint molIndex,
                                     const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;

  for (uint i = 0; i < molKind.nonBonded_1_4.count; ++i) {
    uint p1 = mols.start[molIndex] + molKind.nonBonded_1_4.part1[i];
    uint p2 = mols.start[molIndex] + molKind.nonBonded_1_4.part2[i];
    currentAxes.InRcut(distSq, currentCoords, p1, p2, box);
    if (forcefield.rCutSq > distSq) {
      forcefield.particles->CalcAdd_1_4(energy, distSq,
                                        molKind.AtomKind
                                        (molKind.nonBonded_1_4.part1[i]),
                                        molKind.AtomKind
                                        (molKind.nonBonded_1_4.part2[i]));
      if (electrostatic) {
        qi_qj_Fact = num::qqFact *
                     molKind.AtomCharge(molKind.nonBonded_1_4.part1[i]) *
                     molKind.AtomCharge(molKind.nonBonded_1_4.part2[i]);

        forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
            qi_qj_Fact, false);
      }
    }
  }
}

// Calculate 1-4 nonbonded intra energy using pos
void CalculateEnergy::MolNonbond_1_4(double & energy,
                                    cbmc::TrialMol const &mol,
                                    MoleculeKind const& molKind) const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;
  uint count = molKind.nonBonded_1_4.count;

  for (uint i = 0; i < count; ++i) {
    uint p1 = molKind.nonBonded_1_4.part1[i];
    uint p2 = molKind.nonBonded_1_4.part2[i];
    if(mol.AtomExists(p1) && mol.AtomExists(p2)) {
      currentAxes.InRcut(distSq, mol.GetCoords(), p1, p2, mol.GetBox());
      if (forcefield.rCutSq > distSq) {
        forcefield.particles->CalcAdd_1_4(energy, distSq,
                                          molKind.AtomKind(p1),
                                          molKind.AtomKind(p2));
        if (electrostatic) {
          qi_qj_Fact = num::qqFact * molKind.AtomCharge(p1) *
                      molKind.AtomCharge(p2);

          forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
              qi_qj_Fact, false);
        }
      }
    }
  }
}

// Calculate 1-3 nonbonded intra energy
void CalculateEnergy::MolNonbond_1_3(double & energy,
                                     MoleculeKind const& molKind,
                                     const uint molIndex,
                                     const uint box) const
{
  if (box >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;

  for (uint i = 0; i < molKind.nonBonded_1_3.count; ++i) {
    uint p1 = mols.start[molIndex] + molKind.nonBonded_1_3.part1[i];
    uint p2 = mols.start[molIndex] + molKind.nonBonded_1_3.part2[i];
    currentAxes.InRcut(distSq, currentCoords, p1, p2, box);
    if (forcefield.rCutSq > distSq) {
      forcefield.particles->CalcAdd_1_4(energy, distSq,
                                        molKind.AtomKind
                                        (molKind.nonBonded_1_3.part1[i]),
                                        molKind.AtomKind
                                        (molKind.nonBonded_1_3.part2[i]));
      if (electrostatic) {
        qi_qj_Fact = num::qqFact *
                     molKind.AtomCharge(molKind.nonBonded_1_3.part1[i]) *
                     molKind.AtomCharge(molKind.nonBonded_1_3.part2[i]);

        forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
            qi_qj_Fact, false);
      }
    }
  }
}

// Calculate 1-3 nonbonded intra energy
void CalculateEnergy::MolNonbond_1_3(double & energy, 
                                    cbmc::TrialMol const &mol,
                                    MoleculeKind const& molKind) const
{
  if (mol.GetBox() >= BOXES_WITH_U_B)
    return;

  double distSq;
  double qi_qj_Fact;
  uint count = molKind.nonBonded_1_3.count;

  for (uint i = 0; i < count; ++i) {
    uint p1 = molKind.nonBonded_1_3.part1[i];
    uint p2 = molKind.nonBonded_1_3.part2[i];
    if(mol.AtomExists(p1) && mol.AtomExists(p2)) {
      currentAxes.InRcut(distSq, mol.GetCoords(), p1, p2, mol.GetBox());
      if (forcefield.rCutSq > distSq) {
        forcefield.particles->CalcAdd_1_4(energy, distSq,
                                          molKind.AtomKind(p1),
                                          molKind.AtomKind(p2));
        if (electrostatic) {
          qi_qj_Fact = num::qqFact * molKind.AtomCharge(p1) *
                      molKind.AtomCharge(p2);

          forcefield.particles->CalcCoulombAdd_1_4(energy, distSq,
              qi_qj_Fact, false);
        }
      }
    }
  }
}

// Calculate 1-3 nonbonded intra energy
double CalculateEnergy::IntraEnergy_1_3(const double distSq, const uint atom1,
                                        const uint atom2, const uint molIndex) const
{
  if(!forcefield.OneThree)
    return 0.0;
  else if(forcefield.rCutSq < distSq)
    return 0.0;

  double eng = 0.0;

  MoleculeKind const& thisKind = mols.GetKind(molIndex);
  uint kind1 = thisKind.AtomKind(atom1);
  uint kind2 = thisKind.AtomKind(atom2);

  if (electrostatic) {
    double qi_qj_Fact =  num::qqFact * thisKind.AtomCharge(atom1) *
                         thisKind.AtomCharge(atom2);

    forcefield.particles->CalcCoulombAdd_1_4(eng, distSq, qi_qj_Fact, false);
  }
  forcefield.particles->CalcAdd_1_4(eng, distSq, kind1, kind2);

  if(isnan(eng))
      eng = num::BIGNUM;

  return eng;

}

// Calculate 1-4 nonbonded intra energy
double CalculateEnergy::IntraEnergy_1_4(const double distSq, const uint atom1,
                                        const uint atom2, const uint molIndex) const
{
  if(!forcefield.OneFour)
    return 0.0;
  else if(forcefield.rCutSq < distSq)
    return 0.0;

  double eng = 0.0;


  MoleculeKind const& thisKind = mols.GetKind(molIndex);
  uint kind1 = thisKind.AtomKind(atom1);
  uint kind2 = thisKind.AtomKind(atom2);

  if (electrostatic) {
    double qi_qj_Fact =  num::qqFact * thisKind.AtomCharge(atom1) *
                         thisKind.AtomCharge(atom2);

    forcefield.particles->CalcCoulombAdd_1_4(eng, distSq, qi_qj_Fact, false);
  }
  forcefield.particles->CalcAdd_1_4(eng, distSq, kind1, kind2);

  if(isnan(eng))
      eng = num::BIGNUM;

  return eng;

}

//!Calculates energy and virial tail corrections for the box
void CalculateEnergy::EnergyCorrection(SystemPotential& pot,
                                       BoxDimensions const& boxAxes,
                                       const uint box) const
{
  if (box < BOXES_WITH_U_NB) {
    double en = 0.0;

    for (uint i = 0; i < mols.GetKindsCount(); ++i) {
      uint numI = molLookup.NumKindInBox(i, box);
      for (uint j = 0; j < mols.GetKindsCount(); ++j) {
        uint numJ = molLookup.NumKindInBox(j, box);
        en += mols.pairEnCorrections[i * mols.GetKindsCount() + j] * numI * numJ
              * boxAxes.volInv[box];
      }
    }
    pot.boxEnergy[box].tc = en;
  }
}

//!Calculates energy corrections for the box
double CalculateEnergy::EnergyCorrection(const uint box,
					 const uint *kCount) const
{
  if (box >= BOXES_WITH_U_NB) {
    return 0.0;
  }

  double tc = 0.0;
  for (uint i = 0; i < mols.kindsCount; ++i) {
    for (uint j = 0; j < mols.kindsCount; ++j) {
      tc += mols.pairEnCorrections[i * mols.kindsCount + j] * 
	kCount[i] * kCount[j] * currentAxes.volInv[box];
    }
  }
  return tc;
}

void CalculateEnergy::VirialCorrection(Virial& virial,
                                      BoxDimensions const& boxAxes,
                                      const uint box) const
{
  if (box < BOXES_WITH_U_NB) {
    double vir = 0.0;

    for (uint i = 0; i < mols.GetKindsCount(); ++i) {
      uint numI = molLookup.NumKindInBox(i, box);
      for (uint j = 0; j < mols.GetKindsCount(); ++j) {
        uint numJ = molLookup.NumKindInBox(j, box);
        vir += mols.pairVirCorrections[i * mols.GetKindsCount() + j] *
               numI * numJ * boxAxes.volInv[box];
      }
    }
    virial.tc = vir;
  }
}

//! Calculate Torque
void CalculateEnergy::CalculateTorque(vector<uint>& moleculeIndex,
                                      XYZArray const& coordinates,
                                      XYZArray const& com,
                                      XYZArray const& atomForce,
                                      XYZArray const& atomForceRec,
                                      XYZArray& molTorque,
                                      vector<uint>& moveType,
                                      const uint box)
{
  if(multiParticleEnabled && (box < BOXES_WITH_U_NB)) {
    uint m, p, length, start;
    XYZ tempTorque, distFromCOM;

    // make a pointer to atom force and mol force for openmp
    double *torquex = molTorque.x;
    double *torquey = molTorque.y;
    double *torquez = molTorque.z;
    int torqueCount = molTorque.Count();

    molTorque.Reset();

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(m, p, length, start, distFromCOM, tempTorque) \
reduction(+: torquex[:torqueCount], torquey[:torqueCount], torquez[:torqueCount])
#endif
    for(m = 0; m < moleculeIndex.size(); m++) {
      length = mols.GetKind(moleculeIndex[m]).NumAtoms();
      start = mols.MolStart(moleculeIndex[m]);

      //Only if move is rotation
      if(moveType[moleculeIndex[m]]) {
        // atom iterator
        for(p = start; p < start + length; p++) {
          distFromCOM = coordinates.Difference(p, com, (moleculeIndex[m]));
          distFromCOM = currentAxes.MinImage(distFromCOM, box);
          tempTorque = Cross(distFromCOM, atomForce[p] + atomForceRec[p]);
          
          torquex[moleculeIndex[m]] += tempTorque.x;
          torquey[moleculeIndex[m]] += tempTorque.y;
          torquez[moleculeIndex[m]] += tempTorque.z;
        }
      }
    }

  }
}

void CalculateEnergy::ResetForce(XYZArray& atomForce, XYZArray& molForce,
				 uint box)
{
  if(multiParticleEnabled) {
    uint length, start;

    // molecule iterator
    MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);

    while(thisMol != end) {
      length = mols.GetKind(*thisMol).NumAtoms();
      start = mols.MolStart(*thisMol);

      molForce.Set(*thisMol, 0.0, 0.0, 0.0);
      for(uint p = start; p < start + length; p++) {
        atomForce.Set(p, 0.0, 0.0, 0.0);
      }
      thisMol++;
    }
  }
}

bool CalculateEnergy::FindMolInCavity(std::vector< std::vector<uint> > &mol,
                                      const XYZ& center, const XYZ& cavDim,
                                      const XYZArray& invCav, const uint box,
                                      const uint kind, const uint exRatio)
{
  uint k;
  mol.clear();
  mol.resize(molLookup.GetNumKind());
  double maxLength = cavDim.Max();
    
  if(maxLength <= currentAxes.rCut[box])
  {
    CellList::Neighbors n = cellList.EnumerateLocal(center, box);
    while (!n.Done()) {
      if(currentAxes.InCavity(currentCOM.Get(particleMol[*n]), center, cavDim,
			      invCav, box)) {
	uint molIndex = particleMol[*n];
	//if molecule can be transfer between boxes
	if(!molLookup.IsNoSwap(molIndex)) {
	  k = mols.GetMolKind(molIndex);
	  bool exist = std::find(mol[k].begin(), mol[k].end(), molIndex) !=
	    mol[k].end();
	  if(!exist)
	    mol[k].push_back(molIndex);
	}
      }
      n.Next();
    }
  } else {
    MoleculeLookup::box_iterator n = molLookup.BoxBegin(box);
    MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
    while (n != end) {
      if(currentAxes.InCavity(currentCOM.Get(*n), center, cavDim, invCav, box)){
        uint molIndex = *n;
        //if molecule can be transfer between boxes
        if(!molLookup.IsNoSwap(molIndex)) {
          k = mols.GetMolKind(molIndex);
          bool exist = std::find(mol[k].begin(), mol[k].end(), molIndex) !=
            mol[k].end();
          if(!exist)
            mol[k].push_back(molIndex);
        }
      }
      n++;
    }
  }
    
  //If the is exRate and more molecule kind in cavity, return true.
  if(mol[kind].size() >= exRatio)
    return true;
  else
    return false;
}
