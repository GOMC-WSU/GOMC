/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CALCULATEENERGY_H
#define CALCULATEENERGY_H

#include "../lib/BasicTypes.h"
#include "EnergyTypes.h"

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

namespace cbmc { class TrialMol; }

class CalculateEnergy 
{
   public:
      CalculateEnergy(StaticVals const& stat, System const& sys);

      void Init();

      //! Calculates total energy/virial of all boxes in the system
      SystemPotential SystemTotal() const;

      //! Calculates total energy/virial of a single box in the system
      SystemPotential BoxInter(SystemPotential potential,
                               XYZArray const& coords, 
                               XYZArray const& com,
                               BoxDimensions const& boxAxes,
                               const uint box) const;

      //! Calculates intermolecule energy of all boxes in the system
      //! @param coords Particle coordinates to evaluate for
      //! @param boxAxes Box Dimenions to evaluate in
      //! @return System potential assuming no molecule changes
      SystemPotential SystemInter(SystemPotential potential,
                                  XYZArray const& coords, 
                                  XYZArray const& com,
                                  BoxDimensions const& boxAxes) const;

      //! Calculates intermolecular energy of a molecule were it at molCoords
      //! @param molCoords Molecule coordinates
      //! @param molIndex Index of molecule.
      //! @param box Index of box molecule is in. 
      //! @return
      Intermolecular MoleculeInter(XYZArray const& molCoords,
				   const uint molIndex,
                                   const uint box,
				   XYZ const*const newCOM = NULL) const;

      //! Calculates Nonbonded intra energy for candidate positions
      //! @param trialMol Partially built trial molecule.
      //! @param partIndex Index of particle within the molecule
      //! @param trialPos Contains exactly n potential particle positions
      //! @param energy Return array, must be pre-allocated to size n
      //! @param box Index of box molecule is in
      void ParticleNonbonded(double* energy,
                             const cbmc::TrialMol& trialMol,
                             XYZArray const& trialPos,
                             const uint partIndex,
                             const uint box,
                             const uint trials) const;

      //! Calculates Nonbonded intra energy for candidate positions
      //! @param partIndex Index of the particle within the molecule
      //! @param trialPos Array of candidate positions
      //! @param energy Output Array, at least the size of trialpos
      //! @param molIndex Index of molecule
      //! @param box Index of box molecule is in
      void ParticleInter(double* en,
                         XYZArray const& trialPos,
                         const uint partIndex,
                         const uint molIndex,
                         const uint box,
                         const uint trials) const;

      double MoleculeVirial(const uint molIndex,
                            const uint box) const;


      //! Calculates the change in the TC from adding numChange atoms of a kind
      //! @param box Index of box under consideration
      //! @param kind Kind of particle being added or removed
      //! @param add If removed: false (sign=-1); if added: true (sign=+1)
      Intermolecular MoleculeTailChange(const uint box,
                                        const uint kind,
                                        const bool add) const;

      //Calculates intramolecular energy of a full molecule
      void MoleculeIntra(double & bondEn,
                         double & nonBondEn,
                         const uint molIndex,
                         const uint box) const;

  // private: 

      //Calculates full TC for current system
      void FullTailCorrection(SystemPotential& pot, 
                              BoxDimensions const& boxAxes, 
			      const uint box) const;

      //Calculates bond vectors of a full molecule, stores them in vecs
      void BondVectors(XYZArray & vecs,
                       MoleculeKind const& molKind, 
                       const uint molIndex,
                       const uint box) const;

      //Calculates bond stretch intramolecular energy of a full molecule
      void MolBond(double & energy,
                   MoleculeKind const& molKind,
                   XYZArray const& vecs,
                   const uint box) const;

      //Calculates angular bend intramolecular energy of a full molecule
      void MolAngle(double & energy,
                    MoleculeKind const& molKind,
                    XYZArray const& vecs,
                    const uint box) const;

      //Calculates dihedral torsion intramolecular energy of a full molecule
      void MolDihedral(double & energy,
                       MoleculeKind const& molKind,
                       XYZArray const& vecs,
                       const uint box) const;

      //Calculates Nonbonded intramolecule energy of a full molecule
      void MolNonbond(double & energy,
                      MoleculeKind const& molKind,
                      const uint molIndex,
                      const uint box) const;
      
      bool SameMolecule(const uint p1, const uint p2) const
      { return (particleMol[p1] == particleMol[p2]); }

      const Forcefield& forcefield;
      const Molecules& mols;
      const Coordinates& currentCoords;
      const MoleculeLookup& molLookup;
      const BoxDimensions& currentAxes;
      const COM& currentCOM;
      
      std::vector<int> particleKind;
      std::vector<int> particleMol;
};

#endif /*ENERGY_H*/

