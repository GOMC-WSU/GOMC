/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CALCULATE_ENERGY_H
#define CALCULATE_ENERGY_H

#include <vector>

#include "BasicTypes.h"
#include "CellList.h"
#include "EnergyTypes.h"
#include "Ewald.h"
#include "EwaldCached.h"
#include "NoEwald.h"

//
//    CalculateEnergy.h
//    Energy Calculation functions for Monte Carlo simulation
//    Calculates using const references to a particular Simulation's members
//    Brock Jackman Sep. 2013
//
//    Updated to use radial-based intermolecular pressure
//    Jason Mick    Feb. 2014
//
//    Updated with Mohammad B. to use neighbor list.  Fixed various
//    preexisting flaws in neighbor list codebase to make production ready.
//    Jason Mick    Feb. 2015
//

class StaticVals;
class System;
class Forcefield;
class Molecules;
class MoleculeLookup;
class MoleculeKind;
class Coordinates;
class COM;
class XYZArray;
class BoxDimensions;
class Lambda;

namespace cbmc {
class TrialMol;
}

class CalculateEnergy {
public:
  CalculateEnergy(StaticVals &stat, System &sys);

  void UpdateEwald();
  void Init(System &sys);

  //! Calculates total energy/virial of all boxes in the system
  SystemPotential SystemTotal();

  //! Calculates total energy/virial of a single box in the system
  SystemPotential BoxInter(SystemPotential potential, XYZArray const &coords,
                           BoxDimensions const &boxAxes, const uint box);

  //! Calculates force of a single box in the system
  SystemPotential BoxForce(SystemPotential potential, XYZArray const &coords,
                           XYZArray &atomForce, XYZArray &molForce,
                           BoxDimensions const &boxAxes, const uint box);

  //! Calculate force and virial for the box
  Virial VirialCalc(const uint box);

  //! Set the force for atom and mol to zero for box
  void ResetForce(XYZArray &atomForce, XYZArray &molForce, uint box);

  //! Calculates intermolecule energy of all boxes in the system
  //! @param potential Copy of current energy structure to append result to
  //! @param coords Particle coordinates to evaluate for
  //! @param boxAxes Box Dimensions to evaluate in
  //! @return System potential assuming no molecule changes
  SystemPotential SystemInter(SystemPotential potential, XYZArray const &coords,
                              BoxDimensions const &boxAxes);

  //! Calculates intermolecular energy (LJ and coulomb) of a molecule
  //!                           were it at molCoords.
  //! @param potential Copy of current energy structure to append result to
  //! @param molCoords Molecule coordinates
  //! @param molIndex Index of molecule.
  //! @param box Index of box molecule is in.
  //! @param newCOM (optional) If COM has changed for new coordinate,
  //!                          allows for that to be considered.
  bool MoleculeInter(Intermolecular &inter_LJ, Intermolecular &inter_coulomb,
                     XYZArray const &molCoords, const uint molIndex,
                     const uint box) const;

  //! Calculates Nonbonded intra energy (LJ and coulomb )for
  //!                       candidate positions
  //! @param energy Return array, must be pre-allocated to size n
  //! @param trialMol Partially built trial molecule.
  //! @param trialPos Contains exactly n potential particle positions
  //! @param partIndex Index of particle within the molecule
  //! @param box Index of box molecule is in
  //! @param trials Number of trials to loop over in position array. (cbmc)
  void ParticleNonbonded(double *inter, const cbmc::TrialMol &trialMol,
                         XYZArray const &trialPos, const uint partIndex,
                         const uint box, const uint trials) const;

  //! Calculates Nonbonded inter energy (LJ and coulomb)for
  //!                      candidate positions
  //! @param energy Output Array, at least the size of trialpos
  //! @param trialPos Array of candidate positions
  //! @param partIndex Index of the particle within the molecule
  //! @param molIndex Index of molecule
  //! @param box Index of box molecule is in
  //! @param trials Number of trials to loop over in position array. (cbmc)
  void ParticleInter(double *en, double *real, XYZArray const &trialPos,
                     bool *overlap, const uint partIndex, const uint molIndex,
                     const uint box, const uint trials) const;

  //! Calculates change in the ENergyTC from adding numChange atoms of a kind
  //! @param box Index of box under consideration
  //! @param kind Kind of particle being added or removed
  //! @param add If removed: false (sign=-1); if added: true (sign=+1)
  Intermolecular MoleculeTailChange(const uint box, const uint kind,
                                    const bool add) const;

  //! Calculates change in the Virial TC from adding numChange atoms of a kind
  //! @param box Index of box under consideration
  //! @param kind Kind of particle being added or removed
  //! @param add If removed: false (sign=-1); if added: true (sign=+1)
  Intermolecular MoleculeTailVirChange(const uint box, const uint kind,
                                       const bool add) const;

  //! Calculates the change in the TC from changing the lambdaOld -> lambdaNew
  //! @param box Index of box under consideration
  //! @param kind Kind of particle being transformed in lambda
  double MoleculeTailChange(const uint box, const uint kind,
                            const std::vector<uint> &kCount,
                            const double lambdaOld,
                            const double lambdaNew) const;

  //! Calculates intramolecular energy of a full molecule
  void MoleculeIntra(const uint molIndex, const uint box, double *bondEn) const;

  // used in molecule exchange for calculating bonded and intraNonbonded energy
  Energy MoleculeIntra(cbmc::TrialMol const &mol) const;

  //! Calculates Nonbonded 1_3 intramolecule energy of a full molecule
  // for Martini forcefield
  double IntraEnergy_1_3(const double distSq, const uint atom1,
                         const uint atom2, const uint molIndex) const;

  //! Calculates Nonbonded 1_4 intramolecule energy of a full molecule
  // for Martini forcefield
  double IntraEnergy_1_4(const double distSq, const uint atom1,
                         const uint atom2, const uint molIndex) const;
  //! Calculate Torque
  void CalculateTorque(std::vector<uint> &moleculeIndex,
                       XYZArray const &coordinates, XYZArray const &com,
                       XYZArray const &atomForce, XYZArray const &atomForceRec,
                       XYZArray &molTorque, const uint box);

  // Finding the molecule inside cavity using geometric center
  // and store the molecule Index.
  bool FindMolInCavity(std::vector<std::vector<uint>> &mol, const XYZ &center,
                       const XYZ &cavDim, const XYZArray &invCav,
                       const uint box, const uint kind, const uint exRatio);

  //! Calculates energy corrections for the box
  double EnergyCorrection(const uint box, const uint *kCount) const;

  // Calculate inter energy for single molecule in the system.
  void SingleMoleculeInter(Energy &interEnOld, Energy &interEnNew,
                           const double lambdaOldVDW, const double lambdaNewVDW,
                           const double lambdaOldCoulomb,
                           const double lambdaNewCoulomb, const uint molIndex,
                           const uint box) const;

  // Calculate the change in energy due to lambda
  void EnergyChange(Energy *energyDiff, Energy &dUdL_VDW, Energy &dUdL_Coul,
                    const std::vector<double> &lambda_VDW,
                    const std::vector<double> &lambda_Coul, const uint iState,
                    const uint molIndex, const uint box) const;

private:
  //! Calculates full TC energy for one box in current system
  void EnergyCorrection(SystemPotential &pot, BoxDimensions const &boxAxes,
                        const uint box) const;

  //! Calculates full TC virial for one box in current system
  void VirialCorrection(Virial &virial, BoxDimensions const &boxAxes,
                        const uint box) const;

  //! Calculates bond vectors of a full molecule, stores them in vecs
  void BondVectors(XYZArray &vecs, MoleculeKind const &molKind,
                   const uint molIndex, const uint box) const;

  //! Calculates bond vectors using pos, stores them in vecs
  void BondVectors(XYZArray &vecs, cbmc::TrialMol const &mol,
                   std::vector<bool> &bondExist,
                   MoleculeKind const &molKind) const;

  //! Calculates bond stretch intramolecular energy of a full molecule
  void MolBond(double &energy, MoleculeKind const &molKind,
               XYZArray const &vecs, const uint molIndex, const uint box) const;

  //! Calculates bond stretch intramolecular energy of a non-complete molecule
  void MolBond(double &energy, cbmc::TrialMol const &mol, XYZArray const &vecs,
               std::vector<bool> const &bondExist,
               MoleculeKind const &molKind) const;

  //! Calculates angular bend intramolecular energy of a full molecule
  void MolAngle(double &energy, MoleculeKind const &molKind,
                XYZArray const &vecs, const uint box) const;

  //! Calculates angular bend intramolecular energy of a non-complete molecule
  void MolAngle(double &energy, cbmc::TrialMol const &mol, XYZArray const &vecs,
                std::vector<bool> const &bondExist,
                MoleculeKind const &molKind) const;

  //! Calculates dihedral torsion intramolecular energy of a full molecule
  void MolDihedral(double &energy, MoleculeKind const &molKind,
                   XYZArray const &vecs, const uint box) const;

  //! Calculates dihedral torsion intramolecular energy of a non-complete
  //! molecule
  void MolDihedral(double &energy, cbmc::TrialMol const &mol,
                   XYZArray const &vecs, std::vector<bool> const &bondExist,
                   MoleculeKind const &molKind) const;

  //! Calculates Nonbonded 1_N intramolecule energy of a full molecule
  void MolNonbond(double &energy, MoleculeKind const &molKind,
                  const uint molIndex, const uint box) const;

  //! Calculates Nonbonded 1_N intramolecule energy of a non-complete molecule
  void MolNonbond(double &energy, cbmc::TrialMol const &mol,
                  MoleculeKind const &molKind) const;

  //! Calculates Nonbonded 1_4 intramolecule energy of a full molecule
  void MolNonbond_1_4(double &energy, MoleculeKind const &molKind,
                      const uint molIndex, const uint box) const;

  //! Calculates Nonbonded 1_4 intramolecule energy of a non-complete molecule
  void MolNonbond_1_4(double &energy, cbmc::TrialMol const &mol,
                      MoleculeKind const &molKind) const;

  //! Calculates Nonbonded 1_3 intramolecule energy of a full molecule
  // for Martini forcefield
  void MolNonbond_1_3(double &energy, MoleculeKind const &molKind,
                      const uint molIndex, const uint box) const;

  //! Calculates Nonbonded 1_3 intramolecule energy of a non-complete molecule
  // for Martini forcefield
  void MolNonbond_1_3(double &energy, cbmc::TrialMol const &mol,
                      MoleculeKind const &molKind) const;

  // Calculate the change in LRC for each state
  void ChangeLRC(Energy *energyDiff, Energy &dUdL_VDW,
                 const std::vector<double> &lambda_VDW, const uint iState,
                 const uint molIndex, const uint box) const;

  //! For particles in main coordinates array determines if they belong
  //! to same molecule, using internal arrays.
  bool SameMolecule(const uint p1, const uint p2) const {
    // We don't calculate the energy between two atoms of same molecule
    uint pair1 = particleMol[p1];
    uint pair2 = particleMol[p2];
    return (pair1 == pair2);
  }

  double GetLambdaVDW(uint molA, uint molB, uint box) const;
  double GetLambdaCoulomb(uint molA, uint molB, uint box) const;
  uint NumberOfParticlesInsideBox(uint box);

  const Forcefield &forcefield;
  const Molecules &mols;
  const Coordinates &currentCoords;
  const COM &currentCOM;
  // Used by Wolf Calibration.
  Ewald **ptrToCalcEwald;
  const Ewald *calcEwald;
  const Lambda &lambdaRef;
  XYZArray &atomForceRef;
  XYZArray &molForceRef;
  bool multiParticleEnabled;
  bool electrostatic, ewald;

  // stores kindIndex for each global atom idx
  std::vector<int> particleKind;
  // stores molIndex for each global atom idx
  std::vector<int> particleMol;
  // stores atomIndex for each global atom idx
  std::vector<int> particleIndex;
  // stores charge for each global atom idx
  std::vector<double> particleCharge;
  const MoleculeLookup &molLookup;
  const BoxDimensions &currentAxes;
  const CellList &cellList;
};

#endif /*CALCULATE_ENERGY_H*/
