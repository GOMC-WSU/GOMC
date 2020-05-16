/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
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
#ifdef GOMC_CUDA
#include "CalculateEnergyCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
#define NUMBER_OF_NEIGHBOR_CELL 27

//
//    CalculateEnergy.cpp
//    Energy Calculation functions for Monte Carlo simulation
//    Calculates using const references to a particular Simulation's members
//    Brock Jackman Sep. 2013
//
//    Updated to use radial-based intermolecular pressure
//    Jason Mick    Feb. 2014
//

using namespace geom;

CalculateEnergy::CalculateEnergy(StaticVals & stat, System & sys) :
  forcefield(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
  currentCOM(sys.com),
  atomForceRef(sys.atomForceRef),
  molForceRef(sys.molForceRef),
  lambdaRef(sys.lambdaRef),
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
#ifdef GOMC_CUDA
  InitCoordinatesCUDA(forcefield.particles->getCUDAVars(),
                      currentCoords.Count(), maxAtomInMol, currentCOM.Count());
#endif
}

SystemPotential CalculateEnergy::SystemTotal()
{
  SystemPotential pot =
    SystemInter(SystemPotential(), currentCoords, currentCOM, currentAxes);

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
    #pragma omp parallel for default(shared) private(i, bondEnergy) \
    reduction(+:bondEn, nonbondEn, correction)
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
    pot.boxEnergy[b].correction = correction;

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
                                             BoxDimensions const& boxAxes)
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    //calculate LJ interaction and real term of electrostatic interaction
    potential = BoxInter(potential, coords, boxAxes, b);
    //calculate reciprocate term of electrostatic interaction
    potential.boxEnergy[b].recip = calcEwald->BoxReciprocal(b);
  }

  potential.Total();

  return potential;
}


// Calculate the inter energy for Box. Fractional molecule are not allowed in
// this function. Need to implement the GPU function
SystemPotential CalculateEnergy::BoxInter(SystemPotential potential,
                                          XYZArray const& coords,
                                          BoxDimensions const& boxAxes,
                                          const uint box)
{
  //Handles reservoir box case, returning zeroed structure if
  //interactions are off.
  if (box >= BOXES_WITH_U_NB)
    return potential;

  double tempREn = 0.0, tempLJEn = 0.0;
  double distSq, qi_qj_fact, lambdaVDW, lambdaCoulomb;
  int i;
  XYZ virComponents;

  std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
  std::vector< std::vector<int> > neighborList;
  cellList.GetCellListNeighbor(box, currentCoords.Count(),
      cellVector, cellStartIndex, mapParticleToCell);
  neighborList = cellList.GetNeighborList(box);
  int numberOfCells = neighborList.size();
  int atomNumber = currentCoords.Count();
  int currParticle, currCell, nCellIndex, neighborCell, endIndex, nParticleIndex, nParticle;

#ifdef GOMC_CUDA
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

  CallBoxInterGPU(forcefield.particles->getCUDAVars(), cellVector, cellStartIndex,
      neighborList, mapParticleToCell, coords, boxAxes, electrostatic, particleCharge,
      particleKind, particleMol, REn, LJEn, forcefield.sc_coul,
      forcefield.sc_sigma_6, forcefield.sc_alpha,
      forcefield.sc_power, box);
  tempREn = REn;
  tempLJEn = LJEn;
#else
#ifdef _OPENMP
#pragma omp parallel for default(shared) \
  private(currParticle, currCell, nCellIndex, neighborCell, endIndex, \
      nParticleIndex, nParticle, distSq, qi_qj_fact, virComponents)\
  reduction(+:tempREn, tempLJEn)
#endif
  // loop over all particles
  for(currParticle = 0; currParticle < atomNumber; currParticle++) {
    // find the which cell currParticle belong to
    currCell = mapParticleToCell[currParticle];

    // loop over currCell neighboring cells
    for(nCellIndex = 0; nCellIndex < NUMBER_OF_NEIGHBOR_CELL; nCellIndex++) {
      // find the index of neighboring cell
      neighborCell = neighborList[currCell][nCellIndex];

      // find the ending index in neighboring cell
      endIndex = neighborCell != numberOfCells - 1 ? cellStartIndex[neighborCell+1] : atomNumber;
      // loop over particle inside neighboring cell
      for(nParticleIndex = cellStartIndex[neighborCell];
          nParticleIndex < endIndex; nParticleIndex++) {
        nParticle = cellVector[nParticleIndex];

        // avoid same particles and duplicate work
        if(currParticle < nParticle && particleMol[currParticle] != particleMol[nParticle]) {
          if(boxAxes.InRcut(distSq, virComponents, coords, currParticle, nParticle, box)) {
            lambdaVDW = GetLambdaVDW(particleMol[currParticle], particleMol[nParticle], box);
            if (electrostatic) {
              lambdaCoulomb = GetLambdaCoulomb(particleMol[currParticle],
                  particleMol[nParticle], box);
              qi_qj_fact = particleCharge[currParticle] *
                particleCharge[nParticle] * num::qqFact;
              tempREn += forcefield.particles->CalcCoulomb(distSq,
                  particleKind[currParticle], particleKind[nParticle],
                  qi_qj_fact, lambdaCoulomb, box);
            }
            tempLJEn += forcefield.particles->CalcEn(distSq,
                particleKind[currParticle], particleKind[nParticle], lambdaVDW);
          }
        }
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

SystemPotential CalculateEnergy::BoxForce(SystemPotential potential,
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
  double distSq, qi_qj_fact, lambdaVDW, lambdaCoulomb;
  int i;
  XYZ virComponents, forceLJ, forceReal;
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


  std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
  std::vector<std::vector<int> > neighborList;
  cellList.GetCellListNeighbor(box, currentCoords.Count(), cellVector, cellStartIndex, mapParticleToCell);
  neighborList = cellList.GetNeighborList(box);
  int numberOfCells = neighborList.size();
  int atomNumber = currentCoords.Count();
  int currParticle, currCell, nCellIndex, neighborCell, endIndex, nParticleIndex, nParticle;

#ifdef GOMC_CUDA
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

  CallBoxForceGPU(forcefield.particles->getCUDAVars(), cellVector, 
      cellStartIndex, neighborList, neighborList, mapParticleToCell,
      coords, boxAxes, electrostatic, particleCharge,
      particleKind, particleMol, REn, LJEn, 
      aForcex, aForcey, aForcez, mForcex, mForcey, mForcez,
      atomCount, molCount, reset_force, copy_back,
      forcefield.sc_coul,forcefield.sc_sigma_6, forcefield.sc_alpha,
      forcefield.sc_power, box);

  tempREn = REn;
  tempLJEn = LJEn;

#else
#if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
#pragma omp parallel for default(shared) private(currParticle, currCell, nCellIndex, \
    neighborCell, endIndex, nParticleIndex, nParticle, distSq, qi_qj_fact, \
    virComponents, forceReal, forceLJ, lambdaVDW, lambdaCoulomb) \
  reduction(+:tempREn, tempLJEn, aForcex[:atomCount], aForcey[:atomCount], \
      aForcez[:atomCount], mForcex[:molCount], mForcey[:molCount], mForcez[:molCount])
#endif
  for(currParticle = 0; currParticle < atomNumber; currParticle++) {
    currCell = mapParticleToCell[currParticle];

    for(nCellIndex = 0; nCellIndex < NUMBER_OF_NEIGHBOR_CELL; nCellIndex++) {
      neighborCell = neighborList[currCell][nCellIndex];

      endIndex = neighborCell != numberOfCells - 1 ?
        cellStartIndex[neighborCell+1] : atomNumber;
      for(nParticleIndex = cellStartIndex[neighborCell];
          nParticleIndex < endIndex; nParticleIndex++) {
        nParticle = cellVector[nParticleIndex];

        if(currParticle < nParticle && particleMol[currParticle] != particleMol[nParticle]) {
          if(boxAxes.InRcut(distSq, virComponents, coords, currParticle, nParticle, box)) {
            lambdaVDW = GetLambdaVDW(particleMol[currParticle], particleMol[nParticle], box);

            if (electrostatic) {
              lambdaCoulomb = GetLambdaCoulomb(particleMol[currParticle],
                  particleMol[nParticle], box);
              qi_qj_fact = particleCharge[currParticle] * particleCharge[nParticle] *
                num::qqFact;
              tempREn += forcefield.particles->CalcCoulomb(distSq, particleKind[currParticle],
                  particleKind[nParticle], qi_qj_fact, lambdaCoulomb, box);
            }
            tempLJEn += forcefield.particles->CalcEn(distSq, particleKind[currParticle],
                particleKind[nParticle], lambdaVDW);

            // Calculating the force
            if(multiParticleEnabled) {
              if(electrostatic) {
                forceReal = virComponents *
                  forcefield.particles->CalcCoulombVir(distSq, particleKind[currParticle],
                      particleKind[nParticle], qi_qj_fact, lambdaCoulomb, box);
              }
              forceLJ = virComponents *
                forcefield.particles->CalcVir(distSq, particleKind[currParticle],
                    particleKind[nParticle], lambdaVDW);
              aForcex[currParticle] += forceLJ.x + forceReal.x;
              aForcey[currParticle] += forceLJ.y + forceReal.y;
              aForcez[currParticle] += forceLJ.z + forceReal.z;
              aForcex[nParticle] += -(forceLJ.x + forceReal.x);
              aForcey[nParticle] += -(forceLJ.y + forceReal.y);
              aForcez[nParticle] += -(forceLJ.z + forceReal.z);
              mForcex[particleMol[currParticle]] += (forceLJ.x + forceReal.x);
              mForcey[particleMol[currParticle]] += (forceLJ.y + forceReal.y);
              mForcez[particleMol[currParticle]] += (forceLJ.z + forceReal.z);
              mForcex[particleMol[nParticle]] += -(forceLJ.x + forceReal.x);
              mForcey[particleMol[nParticle]] += -(forceLJ.y + forceReal.y);
              mForcez[particleMol[nParticle]] += -(forceLJ.z + forceReal.z);
            }
          }
        }
      }
    }
  }
#endif
  // setting energy and virial of LJ interaction
  potential.boxEnergy[box].inter = tempLJEn;
  // setting energy and virial of coulomb interaction
  potential.boxEnergy[box].real = tempREn;

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

  double distSq, pVF, pRF, qi_qj, lambdaVDW, lambdaCoulomb;
  int i;
  XYZ virC, comC;

  std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
  std::vector<std::vector<int> > neighborList;
  cellList.GetCellListNeighbor(box, currentCoords.Count(), cellVector,
      cellStartIndex, mapParticleToCell);
  neighborList = cellList.GetNeighborList(box);
  int numberOfCells = neighborList.size();
  int atomNumber = currentCoords.Count();
  int currParticle, currCell, nCellIndex, neighborCell, endIndex, nParticleIndex, nParticle;

#ifdef GOMC_CUDA
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

  CallBoxInterForceGPU(forcefield.particles->getCUDAVars(), 
      cellVector, cellStartIndex, neighborList, mapParticleToCell, 
      currentCoords, currentCOM, currentAxes,
      electrostatic, particleCharge, particleKind,
      particleMol, rT11, rT12, rT13, rT22, rT23, rT33,
      vT11, vT12, vT13, vT22, vT23, vT33,
      forcefield.sc_coul,
      forcefield.sc_sigma_6, forcefield.sc_alpha,
      forcefield.sc_power, box);
#else
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, pVF, pRF, qi_qj, \
    virC, comC, lambdaVDW, lambdaCoulomb) reduction(+:vT11, vT12, vT13, vT22, \
      vT23, vT33, rT11, rT12, rT13, rT22, rT23, rT33)
#endif
  for(currParticle = 0; currParticle < atomNumber; currParticle++) {
    currCell = mapParticleToCell[currParticle];

    for(nCellIndex = 0; nCellIndex < NUMBER_OF_NEIGHBOR_CELL; nCellIndex++) {
      neighborCell = neighborList[currCell][nCellIndex];

      endIndex = neighborCell != numberOfCells - 1 ?
        cellStartIndex[neighborCell+1] : atomNumber;
      for(nParticleIndex = cellStartIndex[neighborCell];
          nParticleIndex < endIndex; nParticleIndex++) {
        nParticle = cellVector[nParticleIndex];

        // make sure the pairs are unique and they belong to different molecules
        if(currParticle < nParticle && particleMol[currParticle] != particleMol[nParticle]) {
          if (currentAxes.InRcut(distSq, virC, currentCoords, currParticle,
                nParticle, box)) {
            pVF = 0.0;
            pRF = 0.0;

            //calculate the distance between com of two molecules
            comC = currentCOM.Difference(particleMol[currParticle],
                particleMol[nParticle]);
            //calculate the minimum image between com of two molecules
            comC = currentAxes.MinImage(comC, box);
            lambdaVDW = GetLambdaVDW(particleMol[currParticle], particleMol[nParticle], box);

            if (electrostatic) {
              lambdaCoulomb = GetLambdaCoulomb(particleMol[currParticle],
                  particleMol[nParticle], box);
              qi_qj = particleCharge[currParticle] * particleCharge[nParticle];

              pRF = forcefield.particles->CalcCoulombVir(distSq, particleKind[currParticle],
                  particleKind[nParticle], qi_qj,
                  lambdaCoulomb, box);
              //calculate the top diagonal of pressure tensor
              rT11 += pRF * (virC.x * comC.x);
              //rT12 += pRF * (0.5 * (virC.x * comC.y + virC.y * comC.x));
              //rT13 += pRF * (0.5 * (virC.x * comC.z + virC.z * comC.x));

              rT22 += pRF * (virC.y * comC.y);
              //rT23 += pRF * (0.5 * (virC.y * comC.z + virC.z * comC.y));

              rT33 += pRF * (virC.z * comC.z);
            }

            pVF = forcefield.particles->CalcVir(distSq, particleKind[currParticle],
                particleKind[nParticle], lambdaVDW);
            //calculate the top diagonal of pressure tensor
            vT11 += pVF * (virC.x * comC.x);
            //vT12 += pVF * (0.5 * (virC.x * comC.y + virC.y * comC.x));
            //vT13 += pVF * (0.5 * (virC.x * comC.z + virC.z * comC.x));

            vT22 += pVF * (virC.y * comC.y);
            //vT23 += pVF * (0.5 * (virC.y * comC.z + virC.z * comC.y));

            vT33 += pVF * (virC.z * comC.z);
          }
        }
      }
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

      double qi_qj_fact, distSq, lambdaVDW, lambdaCoulomb;
      int i;
      XYZ virComponents, forceLJ, forceReal;
      std::vector<uint> nIndex;

      //store atom index in neighboring cell
      while (!n.Done()) {
        nIndex.push_back(*n);
        n.Next();
      }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, \
    virComponents, forceLJ, forceReal, lambdaVDW, lambdaCoulomb) \
      reduction(+:tempREn, tempLJEn)
#endif
      for(i = 0; i < nIndex.size(); i++) {
        distSq = 0.0;
        //Subtract old energy
        if (currentAxes.InRcut(distSq, virComponents, currentCoords, atom,
              nIndex[i], box)) {
          lambdaVDW = GetLambdaVDW(molIndex, particleMol[nIndex[i]], box);

          if (electrostatic) {
            lambdaCoulomb = GetLambdaCoulomb(molIndex, particleMol[nIndex[i]],
                box);
            qi_qj_fact = particleCharge[atom] * particleCharge[nIndex[i]] *
              num::qqFact;

            tempREn -= forcefield.particles->CalcCoulomb(distSq, particleKind[atom],
                particleKind[nIndex[i]], qi_qj_fact, lambdaCoulomb, box);
          }

          tempLJEn -= forcefield.particles->CalcEn(distSq, particleKind[atom],
              particleKind[nIndex[i]], lambdaVDW);
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
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, \
    virComponents, forceReal, forceLJ, lambdaVDW, lambdaCoulomb) \
      reduction(+:tempREn, tempLJEn)
#endif
      for(i = 0; i < nIndex.size(); i++) {
        distSq = 0.0;
        if (currentAxes.InRcut(distSq, virComponents, molCoords, p,
              currentCoords, nIndex[i], box)) {
          lambdaVDW = GetLambdaVDW(molIndex, particleMol[nIndex[i]], box);

          if(distSq < forcefield.rCutLowSq) {
            overlap |= true;
          }

          if (electrostatic) {
            lambdaCoulomb = GetLambdaCoulomb(molIndex, particleMol[nIndex[i]],
                box);
            qi_qj_fact = particleCharge[atom] *
              particleCharge[nIndex[i]] * num::qqFact;

            tempREn += forcefield.particles->CalcCoulomb(distSq,
                particleKind[atom], particleKind[nIndex[i]],
                qi_qj_fact, lambdaCoulomb, box);
          }

          tempLJEn += forcefield.particles->CalcEn(distSq,
              particleKind[atom],
              particleKind[nIndex[i]], lambdaVDW);
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
              kind.AtomKind(*partner), 1.0);
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
  double distSq, qi_qj_Fact, tempLJ, tempReal, lambdaVDW, lambdaCoulomb;
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
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_Fact, \
    lambdaVDW, lambdaCoulomb) reduction(+:tempLJ, tempReal)
#endif
    for(i = 0; i < nIndex.size(); i++) {
      distSq = 0.0;
      if(currentAxes.InRcut(distSq, trialPos, t, currentCoords, nIndex[i], box)) {
        lambdaVDW = GetLambdaVDW(molIndex, particleMol[nIndex[i]], box);

        if(distSq < forcefield.rCutLowSq) {
          overlap[t] |= true;
        }
        tempLJ += forcefield.particles->CalcEn(distSq, kindI,
            particleKind[nIndex[i]],
            lambdaVDW);
        if(electrostatic) {
          lambdaCoulomb = GetLambdaCoulomb(molIndex, particleMol[nIndex[i]],
              box);
          qi_qj_Fact = particleCharge[nIndex[i]] * kindICharge * num::qqFact;
          tempReal += forcefield.particles->CalcCoulomb(distSq, kindI,
              particleKind[nIndex[i]], qi_qj_Fact, lambdaCoulomb, box);
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
      double rhoDeltaIJ_2 = sign * 2.0 * (double)(molLookup.NumKindInBox(j, box)) *
        currentAxes.volInv[box];
      delta.energy += mols.pairEnCorrections[mkIdxIJ] * rhoDeltaIJ_2;
    }

    //We already calculated part of the change for this type in the loop
    delta.energy += mols.pairEnCorrections[mkIdxII] *
      currentAxes.volInv[box];
  }
  return delta;
}

//Calculates the change in the Virial TC from adding numChange atoms of a kind
Intermolecular CalculateEnergy::MoleculeTailVirChange(const uint box,
    const uint kind,
    const bool add) const
{
  Intermolecular delta;

  if (box < BOXES_WITH_U_NB) {

    double sign = (add ? 1.0 : -1.0);
    uint mkIdxII = kind * mols.GetKindsCount() + kind;
    for (uint j = 0; j < mols.GetKindsCount(); ++j) {
      uint mkIdxIJ = j * mols.GetKindsCount() + kind;
      double rhoDeltaIJ_2 = sign * 2.0 * (double)(molLookup.NumKindInBox(j, box)) *
        currentAxes.volInv[box];
      delta.virial += mols.pairVirCorrections[mkIdxIJ] * rhoDeltaIJ_2;
    }

    //We already calculated part of the change for this type in the loop
    delta.virial += mols.pairVirCorrections[mkIdxII] *
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
          (molKind.nonBonded.part2[i]), 1.0);
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
            molKind.AtomKind(p2), 1.0);
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

//!Calculates energy tail corrections for the box
void CalculateEnergy::EnergyCorrection(SystemPotential& pot,
    BoxDimensions const& boxAxes,
    const uint box) const
{
  if(box >= BOXES_WITH_U_NB) {
    return;
  }

  double en = 0.0;
  for (uint i = 0; i < mols.GetKindsCount(); ++i) {
    uint numI = molLookup.NumKindInBox(i, box);
    for (uint j = 0; j < mols.GetKindsCount(); ++j) {
      uint numJ = molLookup.NumKindInBox(j, box);
      en += mols.pairEnCorrections[i * mols.GetKindsCount() + j] * numI * numJ
        * boxAxes.volInv[box];
    }
  }

  if(!forcefield.freeEnergy) {
    pot.boxEnergy[box].tc = en;
  }
#if ENSEMBLE == NVT || ENSEMBLE == NPT
  else {
    //Get the kind and lambda value
    uint fk = lambdaRef.GetKind(box);
    double lambdaVDW = lambdaRef.GetLambdaVDW(fk, box);
    //remove the LRC for one molecule with lambda = 1
    en += MoleculeTailChange(box, fk, false).energy;

    //Add the LRC for fractional molecule
    for (uint i = 0; i < mols.GetKindsCount(); ++i) {
      uint molNum = molLookup.NumKindInBox(i, box);
      if(i == fk) {
        --molNum; // We have one less molecule (it is fractional molecule)
      }
      double rhoDeltaIJ_2 = 2.0 * (double)(molNum) * currentAxes.volInv[box];
      en += lambdaVDW * mols.pairEnCorrections[fk * mols.GetKindsCount() + i] *
        rhoDeltaIJ_2;
    }
    //We already calculated part of the change for this type in the loop
    en += lambdaVDW * mols.pairEnCorrections[fk * mols.GetKindsCount() + fk] *
      currentAxes.volInv[box];
    pot.boxEnergy[box].tc = en;
  }
#endif
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
  if(box >= BOXES_WITH_U_NB) {
    return;
  }
  double vir = 0.0;

  for (uint i = 0; i < mols.GetKindsCount(); ++i) {
    uint numI = molLookup.NumKindInBox(i, box);
    for (uint j = 0; j < mols.GetKindsCount(); ++j) {
      uint numJ = molLookup.NumKindInBox(j, box);
      vir += mols.pairVirCorrections[i * mols.GetKindsCount() + j] *
        numI * numJ * boxAxes.volInv[box];
    }
  }

  if(!forcefield.freeEnergy) {
    virial.tc = vir;
  }
#if ENSEMBLE == NVT || ENSEMBLE == NPT
  else {
    //Get the kind and lambda value
    uint fk = lambdaRef.GetKind(box);
    double lambdaVDW = lambdaRef.GetLambdaVDW(fk, box);
    //remove the LRC for one molecule with lambda = 1
    vir += MoleculeTailVirChange(box, fk, false).virial;

    //Add the LRC for fractional molecule
    for (uint i = 0; i < mols.GetKindsCount(); ++i) {
      uint molNum = molLookup.NumKindInBox(i, box);
      if(i == fk) {
        --molNum; // We have one less molecule (it is fractional molecule)
      }
      double rhoDeltaIJ_2 = 2.0 * (double)(molNum) * currentAxes.volInv[box];
      vir += mols.pairVirCorrections[fk * mols.GetKindsCount() + i] * rhoDeltaIJ_2;
    }
    //We already calculated part of the change for this type in the loop
    vir += mols.pairVirCorrections[fk * mols.GetKindsCount() + fk] *
      currentAxes.volInv[box];
    virial.tc = vir;
  }
#endif
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

#if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
#pragma omp parallel for default(shared) private(m, p, length, start, \
    distFromCOM, tempTorque) reduction(+: torquex[:torqueCount], \
      torquey[:torqueCount], torquez[:torqueCount])
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

  if(maxLength <= currentAxes.rCut[box]) {
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
      if(currentAxes.InCavity(currentCOM.Get(*n), center, cavDim, invCav, box)) {
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


void CalculateEnergy::SingleMoleculeInter(Energy &interEnOld,
    Energy &interEnNew,
    const double lambdaOldVDW,
    const double lambdaNewVDW,
    const double lambdaOldCoulomb,
    const double lambdaNewCoulomb,
    const uint molIndex,
    const uint box) const
{
  double tempREnOld = 0.0, tempLJEnOld = 0.0;
  double tempREnNew = 0.0, tempLJEnNew = 0.0;
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
      XYZ virComponents;
      std::vector<uint> nIndex;

      //store atom index in neighboring cell
      while (!n.Done()) {
        if(particleMol[*n] != molIndex) {
          nIndex.push_back(*n);
        }
        n.Next();
      }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_fact, \
    virComponents) reduction(+:tempREnOld, tempLJEnOld, tempREnNew, tempLJEnNew)
#endif
      for(i = 0; i < nIndex.size(); i++) {
        distSq = 0.0;
        if (currentAxes.InRcut(distSq, virComponents, currentCoords, atom,
              nIndex[i], box)) {

          if (electrostatic) {
            qi_qj_fact = particleCharge[atom] * particleCharge[nIndex[i]] *
              num::qqFact;
            tempREnNew += forcefield.particles->CalcCoulomb(distSq, particleKind[atom],
                particleKind[nIndex[i]], qi_qj_fact, lambdaNewCoulomb, box);
            tempREnOld += forcefield.particles->CalcCoulomb(distSq, particleKind[atom],
                particleKind[nIndex[i]], qi_qj_fact, lambdaOldCoulomb, box);
          }

          tempLJEnNew += forcefield.particles->CalcEn(distSq, particleKind[atom],
              particleKind[nIndex[i]], lambdaNewVDW);
          tempLJEnOld += forcefield.particles->CalcEn(distSq, particleKind[atom],
              particleKind[nIndex[i]], lambdaOldVDW);
        }
      }
    }
  }

  interEnNew.inter = tempLJEnNew;
  interEnNew.real = tempREnNew;
  interEnOld.inter = tempLJEnOld;
  interEnOld.real = tempREnOld;
}

double CalculateEnergy::GetLambdaVDW(uint molA, uint molB, uint box) const
{
  double lambda = 1.0;
  lambda *= lambdaRef.GetLambdaVDW(molA, mols.GetMolKind(molA), box);
  lambda *= lambdaRef.GetLambdaVDW(molB, mols.GetMolKind(molB), box);
  return lambda;
}

double CalculateEnergy::GetLambdaCoulomb(uint molA, uint molB, uint box) const
{
  double lambda = 1.0;
  lambda *= lambdaRef.GetLambdaCoulomb(molA, mols.GetMolKind(molA), box);
  lambda *= lambdaRef.GetLambdaCoulomb(molB, mols.GetMolKind(molB), box);
  //no need for sq root for inter energy. Alwayse one of the molecule has
  // lambda 1
  return lambda;
}

//Calculates the change in the TC from adding numChange atoms of a kind
double CalculateEnergy::MoleculeTailChange(const uint box, const uint kind,
    const std::vector <uint> &kCount,
    const double lambdaOld,
    const double lambdaNew) const
{
  if (box >= BOXES_WITH_U_NB) {
    return 0.0;
  }

  double tcDiff = 0.0;
  uint ktot = mols.GetKindsCount();
  for (uint i = 0; i < ktot; ++i) {
    //We should have only one molecule of fractional kind
    double rhoDeltaIJ_2 = 2.0 * (double)(kCount[i]) * currentAxes.volInv[box];
    uint index = kind * ktot + i;
    tcDiff += (lambdaNew - lambdaOld) * mols.pairEnCorrections[index] *
      rhoDeltaIJ_2;
  }
  uint index = kind * ktot + kind;
  tcDiff += (lambdaNew - lambdaOld) * mols.pairEnCorrections[index] *
    currentAxes.volInv[box];

  return tcDiff;
}

//Calculate the change in energy due to lambda
void CalculateEnergy::EnergyChange(Energy *energyDiff, Energy &dUdL_VDW,
    Energy &dUdL_Coul,
    const std::vector<double> &lambda_VDW,
    const std::vector<double> &lambda_Coul,
    const uint iState, const uint molIndex,
    const uint box) const
{
  if (box >= BOXES_WITH_U_NB) {
    return;
  }

  uint length = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  uint lambdaSize = lambda_VDW.size();
  int i, s;
  double energyOldVDW = 0.0, energyOldCoul = 0.0;
  double *tempLJEnDiff = new double [lambdaSize];
  double *tempREnDiff = new double [lambdaSize];
  double dudl_VDW = 0.0, dudl_Coul = 0.0;
  std::fill_n(tempLJEnDiff, lambdaSize, 0.0);
  std::fill_n(tempREnDiff, lambdaSize, 0.0);

  // Calculate the vdw, short range electrostatic energy
  for (uint p = 0; p < length; ++p) {
    uint atom = start + p;
    CellList::Neighbors n = cellList.EnumerateLocal(currentCoords[atom], box);
    n = cellList.EnumerateLocal(currentCoords[atom], box);

    double qi_qj_fact, distSq;
    XYZ virComponents;
    std::vector<uint> nIndex;

    //store atom index in neighboring cell
    while (!n.Done()) {
      if(particleMol[*n] != molIndex) {
        nIndex.push_back(*n);
      }
      n.Next();
    }

#if defined _OPENMP && _OPENMP >= 201511 // check if OpenMP version is 4.5
#pragma omp parallel for default(shared) private(i, s, distSq, qi_qj_fact, \
    virComponents, energyOldVDW, energyOldCoul) reduction(+:dudl_VDW, dudl_Coul, \
      tempREnDiff[:lambdaSize], tempLJEnDiff[:lambdaSize])
#endif
    for(i = 0; i < nIndex.size(); i++) {
      distSq = 0.0;
      if(currentAxes.InRcut(distSq, virComponents, currentCoords, atom,
            nIndex[i], box)) {

        //Calculate the energy of current state
        energyOldVDW = forcefield.particles->CalcEn(distSq, particleKind[atom],
            particleKind[nIndex[i]],
            lambda_VDW[iState]);
        //Calculate du/dl in VDW for current state
        dudl_VDW += forcefield.particles->CalcdEndL(distSq, particleKind[atom],
            particleKind[nIndex[i]], lambda_VDW[iState]);

        if(electrostatic) {
          qi_qj_fact = particleCharge[atom] * particleCharge[nIndex[i]] *
            num::qqFact;
          energyOldCoul = forcefield.particles->CalcCoulomb(distSq, particleKind[atom],
              particleKind[nIndex[i]], qi_qj_fact,
              lambda_Coul[iState], box);
          //Calculate du/dl in Coulomb for current state.
          dudl_Coul += forcefield.particles->CalcCoulombdEndL(distSq, particleKind[atom],
              particleKind[nIndex[i]], qi_qj_fact,
              lambda_Coul[iState], box);
        }

        for(s = 0; s < lambdaSize; s++) {
          //Calculate the energy of other state
          tempLJEnDiff[s] += forcefield.particles->CalcEn(distSq, particleKind[atom],
              particleKind[nIndex[i]], lambda_VDW[s]);
          tempLJEnDiff[s] += (-energyOldVDW);
          if(electrostatic) {
            tempREnDiff[s] += forcefield.particles->CalcCoulomb(distSq, particleKind[atom],
                particleKind[nIndex[i]], qi_qj_fact,
                lambda_Coul[s], box);
            tempREnDiff[s] += (-energyOldCoul);
          }
        }
      }
    }
  }

  dUdL_VDW.inter = dudl_VDW;
  dUdL_Coul.real = dudl_Coul;
  for(s = 0; s < lambdaSize; s++) {
    energyDiff[s].inter += tempLJEnDiff[s];
    energyDiff[s].real += tempREnDiff[s];
  }
  delete [] tempLJEnDiff;
  delete [] tempREnDiff;

  if (forcefield.useLRC) {
    //Need to calculate change in LRC
    ChangeLRC(energyDiff, dUdL_VDW, lambda_VDW, iState, molIndex, box);
  }
  //Need to calculate change in self
  calcEwald->ChangeSelf(energyDiff, dUdL_Coul, lambda_Coul, iState, molIndex,
      box);
  //Need to calculate change in correction
  calcEwald->ChangeCorrection(energyDiff, dUdL_Coul, lambda_Coul, iState,
      molIndex, box);
  //Need to calculate change in Reciprocal
  calcEwald->ChangeRecip(energyDiff, dUdL_Coul, lambda_Coul, iState, molIndex,
      box);
}

//Calculate the change in LRC for each state
void CalculateEnergy::ChangeLRC(Energy *energyDiff, Energy &dUdL_VDW,
    const std::vector<double> &lambda_VDW,
    const uint iState, const uint molIndex,
    const uint box) const
{
  //Get the kind and lambda value
  uint fk = mols.GetMolKind(molIndex);
  double lambda_istate = lambda_VDW[iState];

  //Add the LRC for fractional molecule
  for(uint s = 0; s < lambda_VDW.size(); s++) {
    double lambdaVDW = lambda_VDW[s];
    for (uint i = 0; i < mols.GetKindsCount(); ++i) {
      uint molNum = molLookup.NumKindInBox(i, box);
      if(i == fk) {
        --molNum; // We have one less molecule (it is fractional molecule)
      }
      double rhoDeltaIJ_2 = 2.0 * (double)(molNum) * currentAxes.volInv[box];
      energyDiff[s].tc += mols.pairEnCorrections[fk * mols.GetKindsCount() + i] *
        rhoDeltaIJ_2 * (lambdaVDW - lambda_istate);
      if(s == iState) {
        //Calculate du/dl in VDW LRC for current state
        dUdL_VDW.tc += mols.pairEnCorrections[fk * mols.GetKindsCount() + i] *
          rhoDeltaIJ_2;
      }
    }
    energyDiff[s].tc += mols.pairEnCorrections[fk * mols.GetKindsCount() + fk] *
      currentAxes.volInv[box] * (lambdaVDW - lambda_istate);
    if(s == iState) {
      //Calculate du/dl in VDW LRC for current state
      dUdL_VDW.tc += mols.pairEnCorrections[fk * mols.GetKindsCount() + fk] *
        currentAxes.volInv[box];
    }
  }
}
