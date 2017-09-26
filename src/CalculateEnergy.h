/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CALCULATEENERGY_H
#define CALCULATEENERGY_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "EwaldCached.h"
#include "Ewald.h"
#include "NoEwald.h"
#include "CellList.h"

#include <vector>

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

namespace cbmc
{
class TrialMol;
}

class CalculateEnergy
{
public:
  CalculateEnergy(StaticVals & stat, System & sys);

  void Init(System & sys);

  //! Calculates total energy/virial of all boxes in the system
  SystemPotential SystemTotal() ;

  //! Calculates total energy/virial of a single box in the system
  SystemPotential BoxInter(SystemPotential potential,
                           XYZArray const& coords,
                           XYZArray const& com,
                           BoxDimensions const& boxAxes,
                           const uint box) ;

  //! Calculate force and virial for the box
  Virial ForceCalc(const uint box);


  //! Calculates intermolecule energy of all boxes in the system
  //! @param potential Copy of current energy structure to append result to
  //! @param coords Particle coordinates to evaluate for
  //! @param com Molecular centers of mass of system under evaluation
  //! @param boxAxes Box Dimenions to evaluate in
  //! @return System potential assuming no molecule changes
  SystemPotential SystemInter(SystemPotential potential,
                              XYZArray const& coords,
                              XYZArray const& com,
                              BoxDimensions const& boxAxes) ;

  //! Calculates intermolecular energy (LJ and coulomb) of a molecule
  //!                           were it at molCoords.
  //! @param potential Copy of current energy structure to append result to
  //! @param molCoords Molecule coordinates
  //! @param molIndex Index of molecule.
  //! @param box Index of box molecule is in.
  //! @param newCOM (optional) If COM has changed for new coordinate,
  //!                          allows for that to be considered.
  void MoleculeInter(Intermolecular &inter_LJ, Intermolecular &inter_coulomb,
                     XYZArray const& molCoords, const uint molIndex,
                     const uint box) const;

  //! Calculates Nonbonded intra energy (LJ and coulomb )for
  //!                       candidate positions
  //! @param energy Return array, must be pre-allocated to size n
  //! @param trialMol Partially built trial molecule.
  //! @param trialPos Contains exactly n potential particle positions
  //! @param partIndex Index of particle within the molecule
  //! @param box Index of box molecule is in
  //! @param trials Number of trials ot loop over in position array. (cbmc)
  void ParticleNonbonded(double* inter, const cbmc::TrialMol& trialMol,
                         XYZArray const& trialPos,
                         const uint partIndex,
                         const uint box,
                         const uint trials) const;


  //! Calculates Nonbonded inter energy (LJ and coulomb)for
  //!                      candidate positions
  //! @param energy Output Array, at least the size of trialpos
  //! @param trialPos Array of candidate positions
  //! @param partIndex Index of the particle within the molecule
  //! @param molIndex Index of molecule
  //! @param box Index of box molecule is in
  //! @param trials Number of trials ot loop over in position array. (cbmc)
  void ParticleInter(double* en, double *real,
                     XYZArray const& trialPos,
                     const uint partIndex,
                     const uint molIndex,
                     const uint box,
                     const uint trials) const;


  //! Calculates the change in the TC from adding numChange atoms of a kind
  //! @param box Index of box under consideration
  //! @param kind Kind of particle being added or removed
  //! @param add If removed: false (sign=-1); if added: true (sign=+1)
  Intermolecular MoleculeTailChange(const uint box,
                                    const uint kind,
                                    const bool add) const;

  //! Calculates intramolecular energy of a full molecule
  double* MoleculeIntra(const uint molIndex, const uint box) const;

  //! Calculates Nonbonded 1_3 intramolecule energy of a full molecule
  //for Martini forcefield
  double IntraEnergy_1_3(const double distSq, const uint atom1,
                         const uint atom2, const uint molIndex) const;

  //! Calculates Nonbonded 1_4 intramolecule energy of a full molecule
  //for Martini forcefield
  double IntraEnergy_1_4(const double distSq, const uint atom1,
                         const uint atom2, const uint molIndex) const;

private:

  //! Calculates full TC energy for one box in current system
  void EnergyCorrection(SystemPotential& pot, BoxDimensions const& boxAxes,
			const uint box) const;

  //! Calculates full TC virial for one box in current system
  void ForceCorrection(Virial& virial, BoxDimensions const& boxAxes,
		       const uint box) const;


  //! Calculates bond vectors of a full molecule, stores them in vecs
  void BondVectors(XYZArray & vecs,
                   MoleculeKind const& molKind,
                   const uint molIndex,
                   const uint box) const;

  //! Calculates bond stretch intramolecular energy of a full molecule
  void MolBond(double & energy,
               MoleculeKind const& molKind,
               XYZArray const& vecs,
               const uint box) const;

  //! Calculates angular bend intramolecular energy of a full molecule
  void MolAngle(double & energy,
                MoleculeKind const& molKind,
                XYZArray const& vecs,
                const uint box) const;

  //! Calculates dihedral torsion intramolecular energy of a full molecule
  void MolDihedral(double & energy,
                   MoleculeKind const& molKind,
                   XYZArray const& vecs,
                   const uint box) const;

  //! Calculates Nonbonded 1_N intramolecule energy of a full molecule
  void MolNonbond(double & energy,
                  MoleculeKind const& molKind,
                  const uint molIndex,
                  const uint box) const;

  //! Calculates Nonbonded 1_4 intramolecule energy of a full molecule
  void MolNonbond_1_4(double & energy,
                      MoleculeKind const& molKind,
                      const uint molIndex,
                      const uint box) const;

  //! Calculates Nonbonded 1_3 intramolecule energy of a full molecule
  //for Martini forcefield
  void MolNonbond_1_3(double & energy,
                      MoleculeKind const& molKind,
                      const uint molIndex,
                      const uint box) const;

  //! For particles in main coordinates array determines if they belong
  //! to same molecule, using internal arrays.
  bool SameMolecule(const uint p1, const uint p2) const
  {
    //We dont calculate the energy between two atom of same molecule or
    uint pair1 = particleMol[p1];
    uint pair2 = particleMol[p2];
    return (pair1 == pair2);
  }

  const Forcefield& forcefield;
  const Molecules& mols;
  const Coordinates& currentCoords;
  const MoleculeLookup& molLookup;
  const BoxDimensions& currentAxes;
  const COM& currentCOM;
  const EwaldCached  *calcEwald;
  bool electrostatic, ewald;

  std::vector<int> particleKind;
  std::vector<int> particleMol;
  std::vector<double> particleCharge;
  const CellList& cellList;
};

#endif /*ENERGY_H*/
