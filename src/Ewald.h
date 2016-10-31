#ifndef EWALD_H
#define EWALD_H

#include "../lib/BasicTypes.h"
#include "EnergyTypes.h"
#include "EwaldCached.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#include <omp.h>
//
//    Calculating Electrostatic calculation without caching Fourier terms.
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald    
//
//    Developed by Y. Li and Mohammad S. Barhaghi
// 
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
class CalculateEnergy;


class Ewald : public EwaldCached 
{
  //friend class CalculateEnergy;
   public:

   Ewald(StaticVals const& stat, System & sys);

   virtual void Init();

   virtual void AllocMem();

   //initiliazie term used for ewald calculation
   virtual void RecipInit(uint box, BoxDimensions const& boxAxes);
   
   //calculate self term for a box
   virtual double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

   //setup reciprocate term for a box
   virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

   //calculate reciprocate energy term for a box
   virtual double BoxReciprocal(uint box) const;

   //calculate correction term for a molecule
   virtual double MolCorrection(uint molIndex, BoxDimensions const& boxAxes,
				uint box)const;

   //calculate reciprocate term for displacement and rotation move
   virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
				const uint box, XYZ const*const newCOM = NULL);	

   //calculate self term for CBMC algorithm
   virtual void SwapSelf(double *self, uint molIndex, uint partIndex, int box, 
			 uint trials) const;
   
   //calculate correction term for linear molecule CBMC algorithm
   virtual void SwapCorrection(double* energy, const cbmc::TrialMol& trialMol, 
			       XYZArray const& trialPos, const uint partIndex, 
			       const uint box, const uint trials) const; 

   //calculate correction term for branched molecule CBMC algorithm
   virtual void SwapCorrection(double* energy, const cbmc::TrialMol& trialMol,
			       XYZArray *trialPos, const int pickedAtom, 
			       uint *partIndexArray, const uint box,
			       const uint trials, const uint PrevIndex,
			       bool Prev) const;  

   //calculate correction term for old configuration
   virtual double CorrectionOldMol(const cbmc::TrialMol& oldMol,
				   const double distSq,
				   const uint i, const uint j) const;

   //calculate reciprocate term in destination box for swap move
   virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box, 
				const int sourceBox, const int molIndex);	

   //calculate reciprocate term in source box for swap move
   virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
				  const uint box, const int molIndex);


   //back up reciptocate value to Ref (will be called during initialization)
   virtual void SetRecipRef(uint box);

   //update reciprocate values
   virtual void UpdateRecip(uint box);

   //restore cosMol and sinMol
   virtual void RestoreMol(int molIndex);

   //update sinMol and cosMol
   virtual void exgMolCache();

};



#endif /*EWALD_H*/
