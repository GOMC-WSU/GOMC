#ifndef CALCULATEENERGY_H
#define CALCULATEENERGY_H

#define CELL_LIST

#include "../lib/BasicTypes.h"
#include "EnergyTypes.h"
#include "cbmc/TrialMol.h"
#ifdef CELL_LIST
#include "CellList.h"
#endif
#include <cmath>
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

namespace cbmc { class TrialMol; }
namespace config_setup{ class SystemVals; }

class CalculateEnergy 
{
   public:
      CalculateEnergy(StaticVals const& stat, System & sys);
      ~CalculateEnergy();

      void Init(config_setup::SystemVals const& val);

      //! Calculates total energy/virial of all boxes in the system
      SystemPotential SystemTotal();
      SystemPotential SystemTotalRecalc(SystemPotential potential,
					XYZArray const& coords,
					XYZArray const& com,
					BoxDimensions const& boxAxes);

      SystemPotential BoxTotalRecalc(SystemPotential potential,
				     XYZArray const& coords,
				     XYZArray const& com,
				     BoxDimensions const& boxAxes,
				     const uint box);

      //! Calculates total energy/virial of a single box in the system
      SystemPotential BoxInter(SystemPotential potential,
                               XYZArray const& coords, 
                               XYZArray const& com,
                               BoxDimensions const& boxAxes,
                               const uint box);

      //! Calculates intermolecule energy of all boxes in the system
      //! @param potential Copy of current energy structure to append result to
      //! @param coords Particle coordinates to evaluate for
      //! @param com Molecular centers of mass of system under evaluation
      //! @param boxAxes Box Dimenions to evaluate in
      //! @return System potential assuming no molecule changes
      SystemPotential SystemInter(SystemPotential potential,
                                  XYZArray const& coords, 
                                  XYZArray const& com,
                                  BoxDimensions const& boxAxes);

      //! Calculates intermolecular energy of a molecule were it at molCoords
      //! @param potential Copy of current energy structure to append result to
      //! @param molCoords Molecule coordinates
      //! @param molIndex Index of molecule.
      //! @param box Index of box molecule is in. 
      //! @param newCOM (optional) If COM has changed for new coordinate,
      //!                          allows for that to be considered.
      //! @return Intermolecular with interactions for that molecule
      void MoleculeInter(Intermolecular &inter,
						Elect &real,
						XYZArray const& molCoords,
						const uint molIndex,
						const uint box,
						XYZ const*const newCOM = NULL) const;

      //! Calculates Nonbonded intra energy for candidate positions
      //! @param energy Return array, must be pre-allocated to size n
      //! @param trialMol Partially built trial molecule.
      //! @param trialPos Contains exactly n potential particle positions
      //! @param partIndex Index of particle within the molecule
      //! @param box Index of box molecule is in
      //! @param trials Number of trials ot loop over in position array. (cbmc)
      void ParticleNonbonded(double* energy,
                             const cbmc::TrialMol& trialMol,
                             XYZArray const& trialPos,
                             const uint partIndex,
                             const uint box,
                             const uint trials) const;

      //! Calculates Nonbonded intra energy for candidate positions
      //! @param energy Output Array, at least the size of trialpos
      //! @param trialPos Array of candidate positions
      //! @param partIndex Index of the particle within the molecule
      //! @param molIndex Index of molecule
      //! @param box Index of box molecule is in
      //! @param trials Number of trials ot loop over in position array. (cbmc)
      void ParticleInter(double* en,
						 double *real, 
                         XYZArray const& trialPos,
                         const uint partIndex,
                         const uint molIndex,
                         const uint box,
                         const uint trials) const;

      //! For insertion moves we calculate the virial only if we accept, to
      //! save work.
      double MoleculeVirial(const uint molIndex,
                            const uint box) const;


      //! Calculates the change in the TC from adding numChange atoms of a kind
      //! @param box Index of box under consideration
      //! @param kind Kind of particle being added or removed
      //! @param add If removed: false (sign=-1); if added: true (sign=+1)
      Intermolecular MoleculeTailChange(const uint box,
                                        const uint kind,
                                        const bool add) const;

      //! Calculates intramolecular energy of a full molecule
      void MoleculeIntra(double & bondEn,
                         double & nonBondEn,
                         const uint molIndex,
                         const uint box) const;

	  //! Ewald starts from here
	  //ewald inter-molecule calculations (real) are done within intermolecular functions

	  //ewald self energy and correction term functions
	  void MolCorrection(double &correction, uint thisMol, int box)const;
	  void BoxSelf(double& self, int box);
	  void SwapSelf(double *self, uint molIndex, uint partIndex, int box, uint trials) const;
	  void SwapCorrection(double* energy, const cbmc::TrialMol& trialMol, XYZArray const& trialPos, const uint partIndex, const uint box, const uint trials) const;

	  //ewald reciprocal energy functions
	  double BoxReciprocal(int box);
	  double MolReciprocal(XYZArray const& molCoords,
							const uint molIndex,
							const uint box,
							XYZ const*const newCOM = NULL);			//calculates difference between new recip and old recip
	  double SwapDestRecip(cbmc::TrialMol Mol, const uint box);		//calculates new recip or old recip only
	  double SwapSourceRecip(uint molIndex, const uint box);

	  void UpdateRestoreRecip(){
	    double *tS, *tC;
	    tS = RecipSinSum;
	    RecipSinSum = SinSumNew;
	    SinSumNew = tS;
	    tC = RecipCosSum;
	    RecipCosSum = CosSumNew;
	    CosSumNew = tC;
	  };

	  //ewald system 

	  //error function
	  //	  double erf(double erfc)const;

	  //const parameter set and retrive functions
	  void Calp(int box, double boxSize) {calp[box] = alpha / boxSize;};
	  void CalpBackup(int box){calp_old[box] = calp[box];};
	  void CalpRestore(){
	    double *swap;
	    swap = calp;
	    calp = calp_old;
	    calp_old = swap;
	  };
	  double kxyz[2][257][3];				//recip vector
	  double prefact[2][257];			//recip variable
	  void SetupRecip(int box);		//set up recip size, recip varaibles, recip vector, and the cache arrays of recip's sin sum and cos sum
	  bool ifEwald(){return DoEwald;};
	  
   private: 

      //! Calculates full TC for one box in current system
      void FullTailCorrection(SystemPotential& pot, 
                              BoxDimensions const& boxAxes, 
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

      //! Calculates Nonbonded intramolecule energy of a full molecule
      void MolNonbond(double & energy,
                      MoleculeKind const& molKind,
                      const uint molIndex,
                      const uint box) const;
  
      //! For particles in main coordinates array determines if they belong
      //! to same molecule, using internal arrays.
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
#ifdef CELL_LIST
      const CellList& cellList;
#endif

      //for ewald sum
      bool DoEwald;
      double kmax1;
      double alpha;
      double *calp, *calp_old;                //alpha over box size
      double *MolSelfEnergy;      //cache self energy for each molecule kind
      int RecipSize[2];
      SystemPotential &sysPotRef;
      std::vector<double> particleCharge;
      //for ewald summation
      double *RecipSinSum;			//RecipSinSum[box][RecipSize]
      double *RecipCosSum;
      double *SinSumNew;
      double *CosSumNew;
};

#endif /*ENERGY_H*/
