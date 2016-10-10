#include "CalculateEnergy.h"        //header for this
#include "EwaldCached.h"            //for ewald calculation
#include "Ewald.h"                  //for ewald calculation
#include "NoEwald.h"                //for ewald calculation
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "../lib/BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"             //
#include "MoleculeLookup.h"
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "cbmc/TrialMol.h"
#include "../lib/GeomLib.h"
#include "../lib/NumLib.h"
#include <cassert>
#include <omp.h>

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

CalculateEnergy::CalculateEnergy(StaticVals const& stat, System & sys) :
   forcefield(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
   currentCOM(sys.com),
#ifdef VARIABLE_PARTICLE_NUMBER
   molLookup(sys.molLookup),
#else
   molLookup(stat.molLookup),
#endif
#ifdef VARIABLE_VOLUME
   currentAxes(sys.boxDimensions)
#else
   currentAxes(stat.boxDimensions)
#endif
   , cellList(sys.cellList) {}


void CalculateEnergy::Init(System & sys)
{
   calcEwald = sys.GetEwald();
   electrostatic = forcefield.electrostatic;
   ewald = forcefield.ewald;
   for(uint m = 0; m < mols.count; ++m)
   {
      const MoleculeKind& molKind = mols.GetKind(m);
      for(uint a = 0; a < molKind.NumAtoms(); ++a)
      {
         particleKind.push_back(molKind.AtomKind(a));
         particleMol.push_back(m);
	 particleCharge.push_back(molKind.AtomCharge(a));
      }
   }
}

SystemPotential CalculateEnergy::SystemTotal() 
{
   SystemPotential pot =
     SystemInter(SystemPotential(), currentCoords, currentCOM, currentAxes);
   //system intra
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      double bondEn = 0.0, nonbondEn = 0.0, self = 0.0, correction = 0.0;
      MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(b),
	end = molLookup.BoxEnd(b);

      while (thisMol != end)
      {
         MoleculeIntra(bondEn, nonbondEn, *thisMol, b);
	 //calculate correction term of electrostatic interaction
	 correction += calcEwald->MolCorrection(*thisMol, currentAxes, b);
	
         ++thisMol;
      }

      pot.boxEnergy[b].intraBond = bondEn;
      pot.boxEnergy[b].intraNonbond = nonbondEn; 
      //calculate self term of electrostatic interaction
      pot.boxEnergy[b].self = calcEwald->BoxSelf(currentAxes, b);
      pot.boxEnergy[b].correction = -1 * correction * num::qqFact; 
      
   }
   
   pot.Total();
   return pot;
}


SystemPotential CalculateEnergy::SystemInter
(SystemPotential potential,
 XYZArray const& coords,
 XYZArray const& com,
 BoxDimensions const& boxAxes) 
{

   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      //calculate LJ interaction and real term of electrostatic interaction
      potential = BoxInter(potential, coords, com, boxAxes, b);
      //calculate reciprocate term of electrostatic interaction
      potential.boxEnergy[b].recip = calcEwald->BoxReciprocal(b);
   }

   potential.Total();
   return potential;
}



SystemPotential CalculateEnergy::BoxInter(SystemPotential potential,
                                          XYZArray const& coords,
                                          XYZArray const& com,
                                          BoxDimensions const& boxAxes,
                                          const uint box)
{
   //Handles reservoir box case, returning zeroed structure if
   //interactions are off.
   if (box >= BOXES_WITH_U_NB) return potential;

   double tempRForce = 0.0, tempLJForce = 0.0, tempREn = 0.0, tempLJEn = 0.0;
   double boxSize = boxAxes.axis.BoxSize(box);
   double distSq, partVirial, partRealVirial, qi_qj_fact, rEn, ljEn;
   uint i;
   XYZ virComponents;
   std::vector<uint> pair1, pair2;
   CellList::Pairs pair = cellList.EnumeratePairs(box);
   //store atom pair index
   while (!pair.Done()) 
   {
     pair1.push_back(pair.First());
     pair2.push_back(pair.Second());
     pair.Next();
   }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, partVirial, partRealVirial, qi_qj_fact, rEn, ljEn, virComponents) reduction(+:tempREn, tempLJEn) reduction(-:tempRForce, tempLJForce) 
#endif    
   for (i = 0; i < pair1.size(); i++)
   {
      if (!SameMolecule(pair1[i], pair2[i]) &&
	  boxAxes.InRcut(distSq, virComponents, coords, pair1[i],
			 pair2[i], box)) 
      {
	 partVirial = 0.0;
	 partRealVirial = 0.0;
	 rEn = 0.0;
	 ljEn = 0.0;
	 if (electrostatic)
	 {
	   qi_qj_fact = particleCharge[pair1[i]] *
	     particleCharge[pair2[i]] * num::qqFact;		  	

	   forcefield.particles->CalcCoulombAdd(rEn, partRealVirial, 
						distSq, qi_qj_fact);
	   tempREn += rEn;
	   //Add to our pressure components for coulomb interaction.
	   tempRForce -= geom::Dot(virComponents * partRealVirial,
				    boxAxes.MinImage
				    (com.Difference(particleMol[pair1[i]], 
						  particleMol[pair2[i]]), box));
	 }
	 
	 forcefield.particles->CalcAdd(ljEn, partVirial, distSq,
				       particleKind[pair1[i]],
				       particleKind[pair2[i]]);
	 tempLJEn += ljEn;

	 //Add to our pressure components for LJ interaction.
	 tempLJForce -= geom::Dot(virComponents * partVirial,
				  boxAxes.MinImage
				  (com.Difference(particleMol[pair1[i]], 
						  particleMol[pair2[i]]), box));
	 
	 
      }      
   }

   
   // setting energy and virial of LJ interaction
   potential.boxEnergy[box].inter = tempLJEn;   
   potential.boxVirial[box].inter = tempLJForce;
   // setting energy and virial of coulomb interaction
   potential.boxEnergy[box].real = tempREn;
   potential.boxVirial[box].real = tempRForce;
   
   if (forcefield.useLRC) 
   {
      FullTailCorrection(potential, boxAxes, box);
   }

   potential.Total();

   return potential;
}

void CalculateEnergy::MoleculeInter(Intermolecular &inter_LJ,
				    Intermolecular &inter_coulomb,
				    XYZArray const& molCoords,
				    const uint molIndex,
				    const uint box,
				    XYZ const*const newCOM) const
{   
   double tempRForce = 0.0, tempLJForce = 0.0, tempREn = 0.0, tempLJEn = 0.0;
   if (box < BOXES_WITH_U_NB)
   {
      double boxSize = currentAxes.axis.BoxSize(box);
      bool hasNewCOM = !(newCOM == NULL);
      uint length = mols.GetKind(molIndex).NumAtoms();
      uint start = mols.MolStart(molIndex);
      double rEn, ljEn;

      for (uint p = 0; p < length; ++p) 
      {
	 uint atom = start + p;
	 CellList::Neighbors n = cellList.EnumerateLocal(currentCoords[atom], 
							 box);
	 n = cellList.EnumerateLocal(currentCoords[atom], box);
	 
	 double partVirial, partRealVirial, qi_qj_fact, distSq; 
	 uint i;
	 XYZ virComponents;
	 std::vector<uint> nIndex;

	 //store atom index in neighboring cell
	 while (!n.Done())
	 {
	   nIndex.push_back(*n);
	   n.Next();
	 }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, partVirial, partRealVirial, qi_qj_fact, rEn, ljEn, virComponents) reduction(+:tempREn, tempLJEn) reduction(-:tempRForce, tempLJForce)
#endif      
	 for(i = 0; i < nIndex.size(); i++)
	 {
	    distSq = 0.0;
	    //Subtract old energy
	    if (currentAxes.InRcut(distSq, virComponents, 
				   currentCoords, atom, nIndex[i], box)) 
	    {
	       partVirial = 0.0;	
	       partRealVirial = 0.0;
	       rEn = 0.0;
	       ljEn = 0.0;
	       if (electrostatic)
	       {
		 qi_qj_fact = particleCharge[atom] * particleCharge[nIndex[i]] *
		   num::qqFact;

		 forcefield.particles->CalcCoulombSub(rEn,
						      partRealVirial, distSq,
						      qi_qj_fact);
		 tempREn += rEn;

		 //Add to our pressure components for coulomb interaction.
		  tempRForce -= geom::Dot(virComponents * partRealVirial,
					  currentAxes.MinImage
					  (currentCOM.Difference
					   (molIndex, particleMol[nIndex[i]]),
					   box));
		 
	       }
		 
	       forcefield.particles->CalcSub(ljEn,
					    partVirial, distSq,
					    particleKind[atom],
					    particleKind[nIndex[i]]);
	       
	       tempLJEn += ljEn;
	       //Add to our pressure components for LJ interaction.
	       tempLJForce -= geom::Dot(virComponents * partVirial,
					currentAxes.MinImage
					(currentCOM.Difference
					 (molIndex, particleMol[nIndex[i]]),
					 box));
	       
	    } 
	 }

	 //add new energy
	 n = cellList.EnumerateLocal(molCoords[p], box);
	 //store atom index in neighboring cell
	 nIndex.clear();
	 while (!n.Done())
	 {
	   nIndex.push_back(*n);
	   n.Next();
	 }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, distSq, partVirial, partRealVirial, qi_qj_fact, rEn, ljEn, virComponents) reduction(+:tempREn, tempLJEn) reduction(-:tempRForce, tempLJForce) 
#endif	 
	 for(i = 0; i < nIndex.size(); i++) 
	 {
	    distSq = 0.0;
	    if (currentAxes.InRcut(distSq, virComponents, 
				   molCoords, p, currentCoords, nIndex[i],
				   box)) 
	    {
	       partVirial = 0.0;	 
	       partRealVirial = 0.0;
	       rEn = 0.0;
	       ljEn = 0.0;
	       if (electrostatic)
	       {
		 qi_qj_fact = particleCharge[atom] *
		   particleCharge[nIndex[i]] * num::qqFact;	       

		 forcefield.particles->CalcCoulombAdd(rEn,
						      partRealVirial, distSq,
						      qi_qj_fact);
		 tempREn += rEn;
		 if (hasNewCOM)
		 {
		   //Add to our pressure components.
		   tempRForce -=
		     geom::Dot(virComponents * partRealVirial,
			       currentAxes.MinImage(*newCOM - currentCOM.Get(particleMol[nIndex[i]]),box));
		   
		 }
		 else
		 {
		   tempRForce -= geom::Dot(virComponents * partRealVirial,
					   currentAxes.MinImage
					   (currentCOM.Difference
					    (molIndex, particleMol[nIndex[i]]),
					    box));
		 }
	       }

	       forcefield.particles->CalcAdd(ljEn, partVirial,
					     distSq, particleKind[atom],
					     particleKind[nIndex[i]]);
	       
	       tempLJEn += ljEn;
	       if (hasNewCOM)
	       {
		  //Add to our pressure components.
		  tempLJForce -= geom::Dot(virComponents * partVirial,
					   currentAxes.MinImage
					   (*newCOM - 
					    currentCOM.Get(particleMol[nIndex[i]]), box));
	       }
	       else
	       {
		  tempLJForce -= geom::Dot(virComponents * partVirial,
					   currentAxes.MinImage
					   (currentCOM.Difference
					    (molIndex, particleMol[nIndex[i]]), box));
	       }
	    }
	 }
      }
   }
   
   inter_LJ.energy = tempLJEn;
   inter_LJ.virial = tempLJForce;
   inter_coulomb.energy = tempREn;
   inter_coulomb.virial = tempRForce;
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

   double boxSize = currentAxes.axis.BoxSize(box);
   const MoleculeKind& kind = trialMol.GetKind();
   //loop over all partners of the trial particle
   const uint* partner = kind.sortedNB.Begin(partIndex);
   const uint* end = kind.sortedNB.End(partIndex);
   while (partner != end)
   {
      if (trialMol.AtomExists(*partner))
      {

         for (uint t = 0; t < trials; ++t)
         {
            double distSq;
	    
            if (currentAxes.InRcut(distSq, trialPos, t, trialMol.GetCoords(),
				   *partner, box))
            {
               inter[t] +=
		 forcefield.particles->CalcEn(distSq,
						    kind.AtomKind(partIndex),
						    kind.AtomKind(*partner));
	       if (electrostatic)
	       {
		 double qi_qj_Fact = kind.AtomCharge(partIndex) *
		   kind.AtomCharge(*partner) * num::qqFact;
		 inter[t] +=
		   forcefield.particles->CalcCoulombEn(distSq, qi_qj_Fact);
	       }
            }
         }
      }
      ++partner;
   }
}

// Calculate 1-4 nonbonded intra energy
// Calculate 1-3 nonbonded intra energy for Martini force field
void CalculateEnergy::ParticleNonbonded_1_4(double* inter,
                                        cbmc::TrialMol const& trialMol,
                                        XYZArray const& trialPos,
                                        const uint partIndex,
                                        const uint box,
                                        const uint trials) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   
   double boxSize = currentAxes.axis.BoxSize(box);
   const MoleculeKind& kind = trialMol.GetKind();


   //loop over all partners of the trial particle
   const uint* partner = kind.sortedNB_1_4.Begin(partIndex);
   const uint* end = kind.sortedNB_1_4.End(partIndex);
   while (partner != end)
   {
      if (trialMol.AtomExists(*partner))
      {

         for (uint t = 0; t < trials; ++t)
         {
            double distSq;
	    if (currentAxes.InRcut(distSq, trialPos, t, trialMol.GetCoords(),
				   *partner, box))
	    {

	      forcefield.particles->CalcAdd_1_4(inter[t], distSq,
						kind.AtomKind(partIndex),
						kind.AtomKind(*partner));
	      if (electrostatic)
	      {
		double qi_qj_Fact = kind.AtomCharge(partIndex) *
		  kind.AtomCharge(*partner) * num::qqFact;
		forcefield.particles->CalcCoulombAdd_1_4(inter[t],
							 distSq, qi_qj_Fact);
	      }
	    }
	    
         }
      }
      ++partner;
   }
}



//! Calculates Nonbonded intra energy for candidate positions in trialPos
void CalculateEnergy::ParticleInter(double* en, double *real,
                                    XYZArray const& trialPos,
                                    const uint partIndex,
                                    const uint molIndex,
                                    const uint box,
                                    const uint trials) const
{
   if (box >= BOXES_WITH_U_NB)
      return;

   double boxSize = currentAxes.axis.BoxSize(box);
   double distSq, qi_qj_Fact, tempLJ, tempReal;
   uint i, t;
   MoleculeKind const& thisKind = mols.GetKind(molIndex);
   uint kindI = thisKind.AtomKind(partIndex);
   double kindICharge = thisKind.AtomCharge(partIndex);
   std::vector<uint> nIndex;

   for(t = 0; t < trials; ++t) 
   {
      nIndex.clear();
      tempReal = 0.0;
      tempLJ = 0.0;
      CellList::Neighbors n = cellList.EnumerateLocal(trialPos[t], box);
      while (!n.Done())
      {
	nIndex.push_back(*n);
	n.Next();
      }

#ifdef _OPENMP 
#pragma omp parallel for default(shared) private(i, distSq, qi_qj_Fact) reduction(+:tempLJ, tempReal)  
#endif   
      for(i = 0; i < nIndex.size(); i++)
      {
	 distSq = 0.0;

         if (currentAxes.InRcut(distSq, trialPos, t, currentCoords,
				nIndex[i], box)) 
	 {
            tempLJ += forcefield.particles->CalcEn(distSq, kindI, 
						 particleKind[nIndex[i]]);
	    if (electrostatic)
	    {
	      qi_qj_Fact = particleCharge[nIndex[i]] * kindICharge * num::qqFact;
	      tempReal += forcefield.particles->CalcCoulombEn(distSq, qi_qj_Fact);
	    }
         }
      }
      en[t] += tempLJ;
      real[t] += tempReal;
   }   

   return;
}

void CalculateEnergy::MoleculeVirial(double & virial_LJ, double & virial_real,
				     const uint molIndex, const uint box) const
{
   double vir_LJ = 0.0;
   double vir_real = 0.0;

   if (box < BOXES_WITH_U_NB)
   {
      double boxSize = currentAxes.axis.BoxSize(box);
      const MoleculeKind& thisKind = mols.GetKind(molIndex);
      uint thisStart = mols.MolStart(molIndex);
      uint thisLength = thisKind.NumAtoms();
      double mag_LJ, mag_real, qi_qj_Fact, distSq;
      uint p, i;
      XYZ virComponents;
      std::vector<uint> nIndex;

      for (p = 0; p < thisLength; ++p) 
      {
	 uint atom = thisStart + p;
	 CellList::Neighbors n =
	    cellList.EnumerateLocal(currentCoords[atom], box);
	 //clear and assign atom Index
	 nIndex.clear();
	 while (!n.Done()) 
	 {
	   nIndex.push_back(*n);
	   n.Next();
	 }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, mag_LJ, mag_real, virComponents) reduction(-:vir_real, vir_LJ)
#endif	 
	 for (i = 0; i < nIndex.size(); i++)
	 {
	    distSq = 0.0;

	    if (particleMol[atom] != particleMol[nIndex[i]] &&
		currentAxes.InRcut(distSq, virComponents, currentCoords, 
				   atom, nIndex[i], box))
	    {
	       mag_LJ = forcefield.particles->CalcVir(distSq,
							 particleKind[atom],
							 particleKind[nIndex[i]]);
	       if (electrostatic)
	       {
		 qi_qj_Fact = particleCharge[atom] * particleCharge[nIndex[i]] *
		   num::qqFact;

		 mag_real =
		   forcefield.particles->CalcCoulombVir(distSq, qi_qj_Fact);

		 vir_real -= geom::Dot(virComponents * mag_real,
				       currentAxes.MinImage
				       (currentCOM.Difference(particleMol[atom],
							      particleMol[nIndex[i]]), box));
	       }
	       
	       vir_LJ -= geom::Dot(virComponents * mag_LJ,
				   currentAxes.MinImage
				   (currentCOM.Difference(particleMol[atom],
							  particleMol[nIndex[i]]), box));
	       
	    }

	 }
      }

   }
   
   virial_LJ = vir_LJ;
   virial_real = vir_real;
}


//Calculates the change in the TC from adding numChange atoms of a kind
Intermolecular CalculateEnergy::MoleculeTailChange(const uint box,
                                                   const uint kind,
                                                   const bool add) const
{
   Intermolecular delta;
   
   if (box < BOXES_WITH_U_NB)
   {
   
      double sign = (add ? 1.0 : -1.0);
      uint mkIdxII = kind * mols.kindsCount + kind;
      for (uint j = 0; j < mols.kindsCount; ++j)
      {
         uint mkIdxIJ = j * mols.kindsCount + kind;
         double rhoDeltaIJ_2 = sign * 2.0 * 
            (double)(molLookup.NumKindInBox(j, box)) * currentAxes.volInv[box];
         delta.energy += mols.pairEnCorrections[mkIdxIJ] * rhoDeltaIJ_2;
         delta.virial -= mols.pairVirCorrections[mkIdxIJ] * rhoDeltaIJ_2;
      }
      //We already calculated part of the change for this type in the loop
      delta.energy += mols.pairEnCorrections[mkIdxII] * 
	 currentAxes.volInv[box];
      delta.virial -= mols.pairVirCorrections[mkIdxII] *
         currentAxes.volInv[box];
   }
   return delta;
}


//Calculates intramolecular energy of a full molecule
void CalculateEnergy::MoleculeIntra(double& bondEn,
                                    double& nonBondEn,
                                    const uint molIndex,
                                    const uint box) const
{
   MoleculeKind& molKind = mols.kinds[mols.kIndex[molIndex]];
   // *2 because we'll be storing inverse bond vectors
   XYZArray bondVec(molKind.bondList.count * 2);
   BondVectors(bondVec, molKind, molIndex, box);
#ifdef _OPENMP
#pragma omp parallel sections 
#endif
{
#ifdef _OPENMP
#pragma omp section
#endif
   MolBond(bondEn, molKind, bondVec, box);

#ifdef _OPENMP
#pragma omp section
#endif
   MolAngle(bondEn, molKind, bondVec, box);

#ifdef _OPENMP
#pragma omp section
#endif
   MolDihedral(bondEn, molKind, bondVec, box);

#ifdef _OPENMP
#pragma omp section
#endif
   MolNonbond(nonBondEn, molKind, molIndex, box);

#ifdef _OPENMP
#pragma omp section
#endif
   MolNonbond_1_4(nonBondEn, molKind, molIndex, box);

#ifdef _OPENMP
#pragma omp section
#endif
   MolNonbond_1_3(nonBondEn, molKind, molIndex, box);
 }
}

void CalculateEnergy::BondVectors(XYZArray & vecs,
                                  MoleculeKind const& molKind,
                                  const uint molIndex,
                                  const uint box) const
{
   for (uint i = 0; i < molKind.bondList.count; ++i)
   {
      uint p1 = mols.start[molIndex] + molKind.bondList.part1[i];
      uint p2 = mols.start[molIndex] + molKind.bondList.part2[i];
      XYZ dist = currentCoords.Difference(p2, p1);
      dist = currentAxes.MinImage(dist, box);

      //store inverse vectors at i+count
      vecs.Set(i, dist);
      vecs.Set(i + molKind.bondList.count, -dist.x, -dist.y, -dist.z);
   }
}


void CalculateEnergy::MolBond(double & energy,
                              MoleculeKind const& molKind,
                              XYZArray const& vecs,
                              const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
 
   for (uint b = 0; b < molKind.bondList.count; ++b)
   {
      energy += forcefield.bonds.Calc(molKind.bondList.kinds[b],
				      vecs.Get(b).Length());
   }
}

void CalculateEnergy::MolAngle(double & energy,
                               MoleculeKind const& molKind,
                               XYZArray const& vecs,
                               const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   for (uint a = 0; a < molKind.angles.Count(); ++a)
   {
      //Note: need to reverse the second bond to get angle properly.
      double theta = Theta(vecs.Get(molKind.angles.GetBond(a, 0)),
			   -vecs.Get(molKind.angles.GetBond(a, 1)));
      energy += forcefield.angles->Calc(molKind.angles.GetKind(a), theta);
   }
}

void CalculateEnergy::MolDihedral(double & energy,
                                  MoleculeKind const& molKind,
                                  XYZArray const& vecs,
                                  const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   for (uint d = 0; d < molKind.dihedrals.Count(); ++d)
   {
      double phi = Phi(vecs.Get(molKind.dihedrals.GetBond(d, 0)),
		       vecs.Get(molKind.dihedrals.GetBond(d, 1)),
		       vecs.Get(molKind.dihedrals.GetBond(d, 2)));
      energy += forcefield.dihedrals.Calc(molKind.dihedrals.GetKind(d), phi);
   }
}

void CalculateEnergy::MolNonbond(double & energy,
                                 MoleculeKind const& molKind,
                                 const uint molIndex,
                                 const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   
   double distSq;
   double qi_qj_Fact;
   double boxSize = currentAxes.axis.BoxSize(box);
   double virial; //we will throw this away
   for (uint i = 0; i < molKind.nonBonded.count; ++i)
   {
      uint p1 = mols.start[molIndex] + molKind.nonBonded.part1[i];
      uint p2 = mols.start[molIndex] + molKind.nonBonded.part2[i];
      
      if (currentAxes.InRcut(distSq, currentCoords, p1, p2, box))
      {
         forcefield.particles->CalcAdd(energy, virial, distSq,
                                      molKind.AtomKind
                                      (molKind.nonBonded.part1[i]),
                                      molKind.AtomKind
                                      (molKind.nonBonded.part2[i]));
	 if (electrostatic)
	 {
	   qi_qj_Fact = num::qqFact * molKind.AtomCharge(molKind.nonBonded.part1[i])
	     * molKind.AtomCharge(molKind.nonBonded.part2[i]);
	 
	   forcefield.particles->CalcCoulombAdd(energy, virial,
						distSq, qi_qj_Fact);
	 }
      }
   }
}

void CalculateEnergy::MolNonbond_1_4(double & energy,
                                 MoleculeKind const& molKind,
                                 const uint molIndex,
                                 const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   
   double distSq;
   double qi_qj_Fact;
   double boxSize = currentAxes.axis.BoxSize(box);
   for (uint i = 0; i < molKind.nonBonded_1_4.count; ++i)
   {
      uint p1 = mols.start[molIndex] + molKind.nonBonded_1_4.part1[i];
      uint p2 = mols.start[molIndex] + molKind.nonBonded_1_4.part2[i];
      if (currentAxes.InRcut(distSq, currentCoords, p1, p2, box)) 
      {
	forcefield.particles->CalcAdd_1_4(energy, distSq,
					  molKind.AtomKind
					  (molKind.nonBonded_1_4.part1[i]),
					  molKind.AtomKind
					  (molKind.nonBonded_1_4.part2[i]));
	if (electrostatic)
	{
	  qi_qj_Fact = num::qqFact *
	    molKind.AtomCharge(molKind.nonBonded_1_4.part1[i]) *
	    molKind.AtomCharge(molKind.nonBonded_1_4.part2[i]);

	  forcefield.particles->CalcCoulombAdd_1_4(energy,
						   distSq, qi_qj_Fact);
	}
      }
   }
}

void CalculateEnergy::MolNonbond_1_3(double & energy,
                                 MoleculeKind const& molKind,
                                 const uint molIndex,
                                 const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   
   double distSq;
   double qi_qj_Fact;
   double boxSize = currentAxes.axis.BoxSize(box);
   for (uint i = 0; i < molKind.nonBonded_1_3.count; ++i)
   {
      uint p1 = mols.start[molIndex] + molKind.nonBonded_1_3.part1[i];
      uint p2 = mols.start[molIndex] + molKind.nonBonded_1_3.part2[i];
      if (currentAxes.InRcut(distSq, currentCoords, p1, p2, box)) 
      {
	forcefield.particles->CalcAdd_1_4(energy, distSq,
					  molKind.AtomKind
					  (molKind.nonBonded_1_3.part1[i]),
					  molKind.AtomKind
					  (molKind.nonBonded_1_3.part2[i]));
	if (electrostatic)
	{
	  qi_qj_Fact = num::qqFact *
	    molKind.AtomCharge(molKind.nonBonded_1_3.part1[i]) *
	    molKind.AtomCharge(molKind.nonBonded_1_3.part2[i]);

	  forcefield.particles->CalcCoulombAdd_1_4(energy,
						   distSq, qi_qj_Fact);
	}
      }
   }
}

double CalculateEnergy::IntraEnergy_1_3(const double distSq, const uint atom1,
				    const uint atom2, const uint molIndex) const
{
   if(!forcefield.OneThree) 
     return 0.0;
   else if(forcefield.rCutSq <= distSq)
     return 0.0;

   double eng = 0.0;
   double  boxSize = 0.0; //dummy because 1-3 is only for Martini

   MoleculeKind const& thisKind = mols.GetKind(molIndex);
   uint kind1 = thisKind.AtomKind(atom1);
   uint kind2 = thisKind.AtomKind(atom2);

   if (electrostatic)
   {
     double qi_qj_Fact =  num::qqFact * thisKind.AtomCharge(atom1) *
       thisKind.AtomCharge(atom2);
    
     forcefield.particles->CalcCoulombAdd_1_4(eng, distSq, qi_qj_Fact);
   }

   forcefield.particles->CalcAdd_1_4(eng, distSq, kind1, kind2);

   return eng;
   
}

double CalculateEnergy::IntraEnergy_1_4(const double distSq, const uint atom1,
				    const uint atom2, const uint molIndex) const
{
   if(!forcefield.OneFour) 
     return 0.0;
   else if(forcefield.rCutSq <= distSq)
     return 0.0;

   double eng = 0.0;
   double  boxSize = 0.0; //dummy because 1-3 is only for Martini

   MoleculeKind const& thisKind = mols.GetKind(molIndex);
   uint kind1 = thisKind.AtomKind(atom1);
   uint kind2 = thisKind.AtomKind(atom2);
   
   if (electrostatic)
   {
     double qi_qj_Fact =  num::qqFact * thisKind.AtomCharge(atom1) *
       thisKind.AtomCharge(atom2);
    
     forcefield.particles->CalcCoulombAdd_1_4(eng, distSq, qi_qj_Fact);
   }
   forcefield.particles->CalcAdd_1_4(eng, distSq, kind1, kind2);

   return eng;
   
}

//!Calculates energy and virial tail corrections for the box
void CalculateEnergy::FullTailCorrection(SystemPotential& pot, 
                                         BoxDimensions const& boxAxes, 
                                         const uint box) const
{
   if (box < BOXES_WITH_U_NB)
   {
      double en = 0.0;
      double vir = 0.0;

      for (uint i = 0; i < mols.kindsCount; ++i)
      {
         uint numI = molLookup.NumKindInBox(i, box);
         for (uint j = 0; j < mols.kindsCount; ++j)
         {
            uint numJ = molLookup.NumKindInBox(j, box);
            en += mols.pairEnCorrections[i * mols.kindsCount + j] * numI * numJ
               * boxAxes.volInv[box];
            vir -= mols.pairVirCorrections[i * mols.kindsCount + j] *
               numI * numJ * boxAxes.volInv[box];
         }
      }
      pot.boxEnergy[box].tc = en;
      pot.boxVirial[box].tc = vir;
   }
}

