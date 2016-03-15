#include "CalculateEnergy.h"        //header for this
#include "Ewald.h"                  //for ewald calculation
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
#ifdef CELL_LIST
   , cellList(sys.cellList)
#endif
{
  calcEwald = sys.GetEwald();
}

void CalculateEnergy::Init()
{
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
	 if (ewald)
	 {
	    correction += calcEwald->MolCorrection(*thisMol, currentAxes, b);
	 }
         ++thisMol;
      }
      pot.boxEnergy[b].intraBond = bondEn;
      pot.boxEnergy[b].intraNonbond = nonbondEn;
      
      if (ewald)
      {
	 pot.boxEnergy[b].self = calcEwald->BoxSelf(currentAxes, b);
	 pot.boxEnergy[b].correction = -1 * correction * num::qqFact;
	 calcEwald->UpdateRecip(b); 
      }
   }
   // if there is volume change in the simulation, initialize cache for backup
   //   if (ewald && BOXES_WITH_U_NB == 2)
   //   {
   //     calcEwald->deAllocateMolCache(&calcEwald->cosMolBoxRecip, &calcEwald->sinMolBoxRecip);
   //     calcEwald->initMolCache(&calcEwald->cosMolBoxRecip, &calcEwald->sinMolBoxRecip);
   //   }
   pot.Total();
   return pot;
}


SystemPotential CalculateEnergy::SystemInter
(SystemPotential potential,
 XYZArray const& coords,
 XYZArray const& com,
 BoxDimensions const& boxAxes) 
{
   if (ewald)
   {
     calcEwald->exgMolCache();
   }
   
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      potential = BoxInter(potential, coords, com, boxAxes, b);
    
      if (ewald)
      {
	calcEwald->RecipInit(b, boxAxes);
	potential.boxEnergy[b].recip = calcEwald->BoxReciprocal(b, coords);
      }
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

   Intermolecular inter, real;
   double boxSize = boxAxes.axis.BoxSize(box);
 
#ifdef CELL_LIST
   CellList::Pairs pair = cellList.EnumeratePairs(box);
   while (!pair.Done()) 
   {
      XYZ virComponents;
      double distSq;
      if (!SameMolecule(pair.First(), pair.Second()) &&
	  boxAxes.InRcut(distSq, virComponents,
			 coords, pair.First(), pair.Second(), box)) 
      {
	 double partVirial = 0.0;
	 double partRealVirial = 0.0;
	 if (electrostatic)
	 {
	   double qi_qj_fact = particleCharge[pair.First()] *
	     particleCharge[pair.Second()] * num::qqFact;		  	

	   forcefield.particles->CalcCoulombAdd(real.energy, partRealVirial, 
						distSq, qi_qj_fact);
	   //Add to our pressure components for coulomb interaction.
	   real.virial -= geom::Dot(virComponents * partRealVirial,
				    boxAxes.MinImage
				    (com.Difference(particleMol[pair.First()], 
						    particleMol[pair.Second()]),
				     box));
	 }
	 
	 forcefield.particles->CalcAdd(inter.energy, partVirial, distSq,
				       particleKind[pair.First()],
				       particleKind[pair.Second()]);

	 //Add to our pressure components for LJ interaction.
	 inter.virial -= geom::Dot(virComponents * partVirial,
				   boxAxes.MinImage
				   (com.Difference(particleMol[pair.First()], 
						   particleMol[pair.Second()]), 
				    box));
	 
	 
      }
      pair.Next();
   }
#else
   MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
   MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
   while (thisMol != end)
   {
      uint m1 = *thisMol;
      MoleculeKind const& thisKind = mols.GetKind(m1);
      //evaluate interaction with all molecules after this
      MoleculeLookup::box_iterator otherMol = thisMol;
      ++otherMol;
      while (otherMol != end)
      {
	 uint m2 = *otherMol;
	 XYZ fMol, fMolReal;
	 MoleculeKind const& otherKind = mols.GetKind(m2);
	 for (uint i = 0; i < thisKind.NumAtoms(); ++i)
         {
	    for (uint j = 0; j < otherKind.NumAtoms(); ++j)
	    {
	       XYZ virComponents;
	       double distSq;
	       if (boxAxes.InRcut(distSq, virComponents,
				  coords, mols.start[m1] + i,
				  mols.start[m2] + j, box))
               {

		  double partVirial = 0.0;
		  double partRealVirial = 0.0;
		  double qi_qj_fact = thisKind.AtomCharge(i) *
		    otherKind.AtomCharge(j) * num::qqFact;
		  
		  forcefield.particles->CalcAdd(inter.energy,
						      partVirial, distSq,
						      thisKind.AtomKind(i),
						      otherKind.AtomKind(j));

		  forcefield.particles->CalcCoulombAdd(real.energy,
						      partRealVirial, distSq,
						      thisKind.AtomKind(i),
						       qi_qj_fact);

		  //Add to our pressure components.
		  fMol += (virComponents * partVirial);
		  fMolReal += (virComponents * partRealVirial);
	       }
	    }
	 }
	 
	 //Pressure is wrt center of masses of two molecules.
	 inter.virial -= geom::Dot(fMol,
			     boxAxes.MinImage
			     (com.Difference(m1, m2), box));

	 real.virial -= geom::Dot(fMolReal,
			     boxAxes.MinImage
			     (com.Difference(m1, m2), box));
	 
	 ++otherMol;
      }
      ++thisMol;
   }
#endif
   
   // setting energy and virial of LJ interaction
   potential.boxEnergy[box].inter = inter.energy;   
   potential.boxVirial[box].inter = inter.virial;
   // setting energy and virial of coulomb interaction
   potential.boxEnergy[box].real = real.energy;
   potential.boxVirial[box].real = real.virial;
   
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
  Intermolecular inter, real;
   if (box < BOXES_WITH_U_NB)
   {
      double boxSize = currentAxes.axis.BoxSize(box);
      bool hasNewCOM = !(newCOM == NULL);
#ifdef CELL_LIST
      uint length = mols.GetKind(molIndex).NumAtoms();
      uint start = mols.MolStart(molIndex);
      for (uint p = 0; p < length; ++p) 
      {
	 uint atom = start + p;
	 CellList::Neighbors n = cellList.EnumerateLocal(currentCoords[atom], 
							 box);
	 while (!n.Done()) 
	 {
	    double distSq = 0.0;
	    XYZ virComponents;
	    //Subtract old energy
	    if (currentAxes.InRcut(distSq, virComponents, 
				   currentCoords, atom, *n, box)) 
	    {
	       double partVirial = 0.0;	
	       double partRealVirial = 0.0;
	       if (electrostatic)
	       {
		 double qi_qj_fact = particleCharge[atom] * particleCharge[*n] *
		   num::qqFact;

		 forcefield.particles->CalcCoulombSub(real.energy,
						      partRealVirial, distSq,
						      qi_qj_fact);
		 //Add to our pressure components for coulomb interaction.
		 		 real.virial -= geom::Dot(virComponents * partRealVirial,
					  currentAxes.MinImage
					  (currentCOM.Difference
					   (molIndex, particleMol[*n]), box));
		 
	       }
		 
	       forcefield.particles->CalcSub(inter.energy,
					    partVirial, distSq,
					    particleKind[atom],
					    particleKind[*n]);
	       
	       //Add to our pressure components for LJ interaction.
	       inter.virial -= geom::Dot(virComponents * partVirial,
					 currentAxes.MinImage
					 (currentCOM.Difference
					  (molIndex, particleMol[*n]), box));
	       
	    } 
	    n.Next();
	 }

	 //add new energy
	 n = cellList.EnumerateLocal(molCoords[p], box);
	 while (!n.Done()) 
	 {
	    double distSq = 0.0;
	    XYZ virComponents;
	    if (currentAxes.InRcut(distSq, virComponents, 
				   molCoords, p, currentCoords, *n, box)) 
	    {
	       double partVirial = 0.0;	 
	       double partRealVirial = 0.0;
	       if (electrostatic)
	       {
		 double qi_qj_fact = particleCharge[atom] * particleCharge[*n] *
		   num::qqFact;	       

		 forcefield.particles->CalcCoulombAdd(real.energy,
						      partRealVirial, distSq,
						      qi_qj_fact);
		 if (hasNewCOM)
		 {
		   //Add to our pressure components.
		   		   real.virial -= geom::Dot(virComponents * partRealVirial,
					    currentAxes.MinImage
					    (*newCOM - 
					     currentCOM.Get(particleMol[*n]),
					     box));
		 }
		 else
		 {
		   real.virial -= geom::Dot(virComponents * partRealVirial,
					    currentAxes.MinImage
					    (currentCOM.Difference
					     (molIndex, particleMol[*n]), 
					     box));
		 }
	       }

	       forcefield.particles->CalcAdd(inter.energy, partVirial,
						   distSq, particleKind[atom],
						   particleKind[*n]);
	       
	       if (hasNewCOM)
	       {
		  //Add to our pressure components.
		  inter.virial -= geom::Dot(virComponents * partVirial,
					     currentAxes.MinImage
					     (*newCOM - 
					      currentCOM.Get(particleMol[*n]),
					      box));
	       }
	       else
	       {
		  inter.virial -= geom::Dot(virComponents * partVirial,
					     currentAxes.MinImage
					     (currentCOM.Difference
					      (molIndex, particleMol[*n]), 
					      box));
	       }
	    }
	    n.Next();
	 }
      }
#else
      MoleculeLookup::box_iterator otherMol = molLookup.BoxBegin(box);
      MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
      uint partStartMolI, partLenMolI, partStartMolJ, partLenMolJ, pkI, pkJ,
	 pOldI, pNewI, pJ;
      partStartMolI = partLenMolI = partStartMolJ = partLenMolJ = pkI = pkJ =
	 pOldI = pNewI = pJ = 0;
      mols.GetRangeStartLength(partStartMolI, partLenMolI, molIndex);
      const MoleculeKind& thisKind = mols.GetKind(molIndex);
      //looping over all molecules in box
      while (otherMol != end)
      {
         uint m2 = *otherMol;
         XYZ fMolO, fMolN, fMolOldReal, fMolNewReal;
         //except itself
         if (m2 != molIndex)
         {
            const MoleculeKind& otherKind = mols.GetKind(m2);
            //compare all particle pairs
            for (uint i = 0; i < partLenMolI; ++i)
            {
               pOldI = i + partStartMolI;
               pkI = thisKind.AtomKind(i);
               mols.GetRangeStartLength(partStartMolJ, partLenMolJ, m2);
               for (uint j = 0; j < partLenMolJ; ++j)
               {
                  XYZ virComponents;
                  double distSq;
                  pJ = j + partStartMolJ;
                  pkJ = otherKind.AtomKind(j);
                  //Subtract old energy
                  if (currentAxes.InRcut(distSq, virComponents,
                                         currentCoords, pOldI, pJ, box))
                  {
                     double partVirial = 0.0;
		     double partRealVirial = 0.0;
		     double qi_qj_fact = thisKind.AtomCharge(i) * 
		       otherKind.AtomCharge(j) * num::qqFact;
                     
                     forcefield.particles->CalcSub
                        (inter.energy, partVirial, distSq, pkI, pkJ);

		     forcefield.particles->CalcCoulombSub(real.energy,
							  partRealVirial,
							  distSq, qi_qj_fact);
                     
                     fMolO += virComponents * partVirial;
		     fMolOldReal += virComponents * partRealVirial;
                  }
                  //Add new energy
                  if (currentAxes.InRcut(distSq, virComponents, molCoords, i,
                                         currentCoords, pJ, box))
                  {
                     double partVirial = 0.0;
                     double partRealVirial = 0.0;
                     forcefield.particles->CalcAdd(inter.energy,
							 partVirial,
							 distSq, pkI, pkJ);
		     forcefield.particles->CalcCoulombAdd(real.energy,
							 partRealVirial,
							  distSq, qi_qj_fact);
                     
                     //Add to our pressure components.
                     fMolN += (virComponents * partVirial);
		     fMolNewReal += (virComponents * partRealVirial);
                  }
               }
            }
            //Pressure is wrt center of masses of two molecules.
            inter.virial -= geom::Dot(fMolO,
                                       currentAxes.MinImage
                                       (currentCOM.Difference(molIndex, m2),
                                        box));

	    real.virial -= geom::Dot(fMolOldReal,
                                       currentAxes.MinImage
                                       (currentCOM.Difference(molIndex, m2),
                                        box));
            if (hasNewCOM)
            {
               inter.virial -= geom::Dot(fMolN,
                                          currentAxes.MinImage
                                          (*newCOM - currentCOM.Get(m2), box));
	       real.virial -= geom::Dot(fMolNewReal,
                                          currentAxes.MinImage
                                          (*newCOM - currentCOM.Get(m2), box));
            }
            else
            {
               inter.virial -= geom::Dot(fMolN,
                                          currentAxes.MinImage
                                          (currentCOM.Difference(molIndex, m2),
                                           box));
	       real.virial -= geom::Dot(fMolNewReal,
                                          currentAxes.MinImage
                                          (currentCOM.Difference(molIndex, m2),
                                           box));
            }
         }
         
         ++otherMol;
      }
#endif
   }
   
   inter_LJ = inter;
   inter_coulomb = real;
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
   double distSq;
   MoleculeKind const& thisKind = mols.GetKind(molIndex);
   uint kindI = thisKind.AtomKind(partIndex);
   double kindICharge = thisKind.AtomCharge(partIndex);
#ifdef CELL_LIST
  //uint kind = particleKind[partIndex];
   for(int t = 0; t < trials; ++t) 
   {
      CellList::Neighbors n = cellList.EnumerateLocal(trialPos[t], box);
      while (!n.Done())
      {
         int atom = *n;
         double distSq = 0.0;

         if (currentAxes.InRcut(distSq, trialPos, t, currentCoords, atom, box)) 
	 {
            en[t] += forcefield.particles->CalcEn(distSq, kindI, 
						 particleKind[atom]);
	    if (electrostatic)
	    {
	      double qi_qj_Fact = particleCharge[atom] * kindICharge *
		num::qqFact;
	      real[t] += forcefield.particles->CalcCoulombEn(distSq, qi_qj_Fact);
	    }
         }
         n.Next();
      }
   }   
#else
   MoleculeLookup::box_iterator molInBox = molLookup.BoxBegin(box);
   MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
   //looping over all molecules in box
   while (molInBox != end)
   {
      uint otherMol = *molInBox;
      //except itself
      if (otherMol != molIndex)
      {
         MoleculeKind const& otherKind = mols.GetKind(otherMol);
         uint otherStart = mols.MolStart(otherMol);
         uint otherLength = otherKind.NumAtoms();
         //compare all particle pairs
         for (uint j = 0; j < otherLength; ++j)
         {
            uint kindJ = otherKind.AtomKind(j);
	    double kindJCharge = otherKind.AtomCharge(j);
	    double qi_qj_Fact = kindICharge * kindJCharge * num::qqFact;

            for (uint t = 0; t < trialPos.Count(); ++t)
            {
	       double distSq = 0.0;
               if (currentAxes.InRcut(distSq, trialPos, t, currentCoords,
                  otherStart + j, box))
               {
                  en[t] += forcefield.particles->CalcEn(distSq, kindI,
							      kindJ);
		  real[t] += forcefield.particles->CalcEn(distSq, qi_qj_Fact);
               }
            }
         }
      }
      ++molInBox;
   }
#endif
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
#ifdef CELL_LIST
      for (uint p = 0; p < thisLength; ++p) 
      {
	 uint atom = thisStart + p;
	 CellList::Neighbors n =
	    cellList.EnumerateLocal(currentCoords[atom], box);
	 while (!n.Done()) 
	 {
	    double distSq = 0.0;
	    XYZ virComponents;

	    if (particleMol[atom] != particleMol[*n] &&
		currentAxes.InRcut(distSq, virComponents, currentCoords, 
				   atom, *n, box))
	    {
	       double mag_LJ = forcefield.particles->CalcVir(distSq,
							 particleKind[atom],
							 particleKind[*n]);
	       if (electrostatic)
	       {
		 double qi_qj_Fact = particleCharge[atom] * particleCharge[*n] *
		   num::qqFact;
		 double mag_real =
		   forcefield.particles->CalcCoulombVir(distSq, qi_qj_Fact);
		 vir_real -= geom::Dot(virComponents * mag_real,
				       currentAxes.MinImage
				       (currentCOM.Difference(particleMol[atom],
							      particleMol[*n]), 
					box));
	       }
	       
	       vir_LJ -= geom::Dot(virComponents * mag_LJ,
				   currentAxes.MinImage
				   (currentCOM.Difference(particleMol[atom],
							  particleMol[*n]), 
				    box));
	       
	    }
	    n.Next();
	 }
      }
#else
      MoleculeLookup::box_iterator molInBox = molLookup.BoxBegin(box);
      MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);      
      //looping over all molecules in box
      while (molInBox != end)
      {
         uint otherMol = *molInBox;
         //except itself
         if(otherMol == molIndex)
         {
            ++molInBox;
            continue;
         }
         const MoleculeKind& otherKind = mols.GetKind(otherMol);
         XYZ forceOnMol, forceOnMolReal;
         for (uint i = 0; i < thisLength; ++i)
         {
            uint kindI = thisKind.AtomKind(i);
	    uint kindICharge = thisKind.AtomCharge(i);
            uint otherStart = mols.MolStart(otherMol);
            uint otherLength = otherKind.NumAtoms();
            for (uint j = 0; j < otherLength; ++j)
            {
               XYZ forceVector;
               double distSq = 0.0;
               uint kindJ = otherKind.AtomKind(j);
	       uint kindJCharge = otherKind.AtomCharge(j);
	       double qi_qj_Fact = kindICharge * kindJCharge * num::qqFact;

               if (currentAxes.InRcut(distSq, forceVector, 
                                      currentCoords, thisStart + i,
                                      otherStart + j, box))
                  {
                     //sum forces between all particles in molecule pair
                     double mag_LJ = forcefield.particles->CalcVir(distSq,
                                                               kindI, kindJ);
		     double mag_real =
		       forcefield.particles->CalcCoulombVir(distSq, qi_qj_Fact);
                     forceOnMol += forceVector * mag_LJ;
		     forceOnMolReal += forceVector * mag_real;
                  }
            }
         }
         //correct for center of mass
         vir_LJ -= geom::Dot(forceOnMol, 
                             currentAxes.MinImage(currentCOM.Get(molIndex) - 
                                                  currentCOM.Get(otherMol),
                                                  box));
	 vir_real -= geom::Dot(forceOnMolReal, 
                             currentAxes.MinImage(currentCOM.Get(molIndex) - 
                                                  currentCOM.Get(otherMol),
                                                  box));
         ++molInBox;
      }
#endif
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
   MolBond(bondEn, molKind, bondVec, box);
   MolAngle(bondEn, molKind, bondVec, box);
   MolDihedral(bondEn, molKind, bondVec, box);
   MolNonbond(nonBondEn, molKind, molIndex, box);
   MolNonbond_1_4(nonBondEn, molKind, molIndex, box);
   MolNonbond_1_3(nonBondEn, molKind, molIndex, box);
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
	 qi_qj_Fact = num::qqFact * molKind.AtomCharge(molKind.nonBonded.part1[i])
	   * molKind.AtomCharge(molKind.nonBonded.part2[i]);
	 
	 forcefield.particles->CalcCoulombAdd(energy, virial,
						  distSq, qi_qj_Fact);
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

   forcefield.particles->CalcAdd_1_4(eng, distSq, kind1, kind2);

   return eng;
   
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

