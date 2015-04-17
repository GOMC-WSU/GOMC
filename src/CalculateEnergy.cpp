/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "CalculateEnergy.h"        //header for this
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

CalculateEnergy::CalculateEnergy(StaticVals const& stat, System const& sys) :
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
{ }

void CalculateEnergy::Init()
{
   for(uint m = 0; m < mols.count; ++m)
   {
      const MoleculeKind& molKind = mols.GetKind(m);
      for(uint a = 0; a < molKind.NumAtoms(); ++a)
      {
         particleKind.push_back(molKind.AtomKind(a));
         particleMol.push_back(m);
      }
   }
}

SystemPotential CalculateEnergy::SystemTotal() const
{
   SystemPotential pot =
      SystemInter(SystemPotential(), currentCoords,
                  currentCOM, currentAxes);
   //system intra
   for (uint box = 0; box < BOX_TOTAL; ++box)
   {
      double bondEn = 0.0, nonbondEn = 0.0;
      MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);
      MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
      while (thisMol != end)
      {
         MoleculeIntra(bondEn, nonbondEn, *thisMol, box);
         ++thisMol;
      }
      pot.boxEnergy[box].intraBond = bondEn;
      pot.boxEnergy[box].intraNonbond = nonbondEn;
   }
   pot.Total();
   return pot;
}

SystemPotential CalculateEnergy::SystemInter
(SystemPotential potential,
 XYZArray const& coords,
 XYZArray const& com,
 BoxDimensions const& boxAxes) const
{
   for (uint box = 0; box < BOXES_WITH_U_NB; ++box)
   {
      potential = BoxInter(potential, coords, com, boxAxes, box);
   }
   potential.Total();
   return potential;
}



SystemPotential CalculateEnergy::BoxInter(SystemPotential potential,
                                          XYZArray const& coords,
                                          XYZArray const& com,
                                          BoxDimensions const& boxAxes,
                                          const uint box) const
{
   //Handles reservoir box case, returning zeroed structure if
   //interactions are off.
   if (box >= BOXES_WITH_U_NB) return potential;

   Intermolecular inter;
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
	 XYZ fMol;
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
		  
		  forcefield.particles.CalcAdd
		     (inter.energy, partVirial, distSq,
		      thisKind.AtomKind(i), otherKind.AtomKind(j));

		  //Add to our pressure components.
		  fMol += (virComponents * partVirial);
	       }
	    }
	 }
	 
	 //Pressure is wrt center of masses of two molecules.
	 inter.virial -= geom::Dot(fMol,
			     boxAxes.MinImage
			     (com.Difference(m1, m2), box));
	 
	 ++otherMol;
      }
      ++thisMol;
   }
   
   potential.boxEnergy[box].inter = inter.energy;
   potential.boxVirial[box].inter = inter.virial;
   
   if (forcefield.useLRC) 
   {
      FullTailCorrection(potential, boxAxes, box);
   }

   potential.Total();

   return potential;
}

Intermolecular CalculateEnergy::MoleculeInter(XYZArray const& molCoords,
                                              const uint molIndex,
                                              const uint box,
                                              XYZ const*const newCOM) const
{
   
   Intermolecular result;
   if (box < BOXES_WITH_U_NB)
   {
      bool hasNewCOM = !(newCOM == NULL);
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
         XYZ fMolO, fMolN;
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
                     
                     forcefield.particles.CalcSub
                        (result.energy, partVirial, distSq, pkI, pkJ);
                     
                     fMolO += virComponents * partVirial;
                  }
                  //Add new energy
                  if (currentAxes.InRcut(distSq, virComponents, molCoords, i,
                                         currentCoords, pJ, box))
                  {
                     double partVirial = 0.0;
                     
                     forcefield.particles.CalcAdd(result.energy, partVirial,
                                                  distSq, pkI, pkJ);
                     
                     //Add to our pressure components.
                     fMolN += (virComponents * partVirial);
                  }
               }
            }
            //Pressure is wrt center of masses of two molecules.
            result.virial -= geom::Dot(fMolO,
                                       currentAxes.MinImage
                                       (currentCOM.Difference(molIndex, m2),
                                        box));
            if (hasNewCOM)
            {
               result.virial -= geom::Dot(fMolN,
                                          currentAxes.MinImage
                                          (*newCOM - currentCOM.Get(m2), box));
            }
            else
            {
               result.virial -= geom::Dot(fMolN,
                                          currentAxes.MinImage
                                          (currentCOM.Difference(molIndex, m2),
                                           box));
            }
         }
         
         ++otherMol;
      }
   }
   return result;
}

void CalculateEnergy::ParticleNonbonded(double* energy,
                                        const cbmc::TrialMol& trialMol,
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
   while (partner != end)
   {
      if (trialMol.AtomExists(*partner))
      {
         for (uint i = 0; i < trialPos.Count(); ++i)
         {
            double distSq;
            if (currentAxes.InRcut(distSq, trialPos, i, trialMol.GetCoords(),
               *partner, box))
            {
               energy[i] += forcefield.particles.CalcEn(distSq,
                  kind.AtomKind(partIndex), kind.AtomKind(*partner));
            }
         }
      }
      ++partner;
   }
}

//! Calculates Nonbonded intra energy for candidate positions in trialPos
void CalculateEnergy::ParticleInter(double* en,
                                    XYZArray const& trialPos,
                                    const uint partIndex,
                                    const uint molIndex,
                                    const uint box,
                                    const uint trials) const
{
   if (box >= BOXES_WITH_U_NB)
      return;
   
   double distSq;
   MoleculeLookup::box_iterator molInBox = molLookup.BoxBegin(box);
   MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
   MoleculeKind const& thisKind = mols.GetKind(molIndex);
   uint kindI = thisKind.AtomKind(partIndex);

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
            for (uint i = 0; i < trialPos.Count(); ++i)
            {
               if (currentAxes.InRcut(distSq, trialPos, i, currentCoords,
                  otherStart + j, box))
               {
                  en[i] += forcefield.particles.CalcEn(distSq, kindI, kindJ);
               }

			
            }
         }
      }
      ++molInBox;
   }

 
    return;
}

double CalculateEnergy::MoleculeVirial(const uint molIndex,
                                       const uint box) const
{
   double virial = 0;
   if (box < BOXES_WITH_U_NB)
   {
      MoleculeLookup::box_iterator molInBox = molLookup.BoxBegin(box);
      MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);
      
      const MoleculeKind& thisKind = mols.GetKind(molIndex);
      uint thisStart = mols.MolStart(molIndex);
      uint thisLength = thisKind.NumAtoms();
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
         XYZ forceOnMol;
         for (uint i = 0; i < thisLength; ++i)
         {
            uint kindI = thisKind.AtomKind(i);
            uint otherStart = mols.MolStart(otherMol);
            uint otherLength = otherKind.NumAtoms();
            for (uint j = 0; j < otherLength; ++j)
            {
               XYZ forceVector;
               double distSq;
               uint kindJ = otherKind.AtomKind(j);
               if (currentAxes.InRcut(distSq, forceVector, 
                                      currentCoords, thisStart + i,
                                      otherStart + j, box))
                  {
                     //sum forces between all particles in molecule pair
                     double mag = forcefield.particles.CalcVir(distSq,
                                                               kindI, kindJ);
                     forceOnMol += forceVector * mag;
                  }
            }
         }
         //correct for center of mass
         virial -= geom::Dot(forceOnMol, 
                             currentAxes.MinImage(currentCOM.Get(molIndex) - 
                                                  currentCOM.Get(otherMol),
                                                  box));
         ++molInBox;
      }
      double what = virial;
   }
   return virial;
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
   for (uint i = 0; i < molKind.bondList.count; ++i)
   {
      energy += forcefield.bonds.Calc(molKind.bondList.kinds[i],
				      vecs.Get(i).Length());
   }
}

void CalculateEnergy::MolAngle(double & energy,
                               MoleculeKind const& molKind,
                               XYZArray const& vecs,
                               const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   for (uint i = 0; i < molKind.angles.Count(); ++i)
   {
      double theta = Theta(vecs.Get(molKind.angles.GetBond(i, 0)),
         -vecs.Get(molKind.angles.GetBond(i, 1)));
      energy += forcefield.angles.Calc(molKind.angles.GetKind(i), theta);
   }
}

void CalculateEnergy::MolDihedral(double & energy,
                                  MoleculeKind const& molKind,
                                  XYZArray const& vecs,
                                  const uint box) const
{
   if (box >= BOXES_WITH_U_B)
      return;
   for (uint i = 0; i < molKind.dihedrals.Count(); ++i)
   {
      double phi = Phi(vecs.Get(molKind.dihedrals.GetBond(i, 0)),
         vecs.Get(molKind.dihedrals.GetBond(i, 1)),
         vecs.Get(molKind.dihedrals.GetBond(i, 2)));
      energy += forcefield.dihedrals.Calc(molKind.dihedrals.GetKind(i), phi);
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
   double virial; //we will throw this away
   for (uint i = 0; i < molKind.nonBonded.count; ++i)
   {
      uint p1 = mols.start[molIndex] + molKind.nonBonded.part1[i];
      uint p2 = mols.start[molIndex] + molKind.nonBonded.part2[i];
      if (currentAxes.InRcut(distSq, currentCoords, p1, p2, box))
      {
         forcefield.particles.CalcAdd(energy, virial, distSq,
                                      molKind.AtomKind
                                      (molKind.nonBonded.part1[i]),
                                      molKind.AtomKind
                                      (molKind.nonBonded.part2[i]));
      }
   }
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


