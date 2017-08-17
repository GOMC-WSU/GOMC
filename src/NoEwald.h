/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.0
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef NOEWALD_H
#define NOEWALD_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "EwaldCached.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif
//
//    Called when Ewlad method was not used.
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


class NoEwald : public EwaldCached
{
  //friend class CalculateEnergy;
   public:

   NoEwald(StaticVals const& stat, System & sys);

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

   //calculate reciprocate force term for a box
   virtual Virial ForceReciprocal(Virial& virial, uint box) const;

   //calculate correction force term for a box
   virtual Virial ForceCorrection(Virial& virial, uint box) const;

   //calculate correction term for a molecule
   virtual double MolCorrection(uint molIndex, uint box)const;

   //calculate reciprocate term for displacement and rotation move
   virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
				const uint box, XYZ const*const newCOM = NULL);

   //calculate self term after swap move
   virtual double SwapSelf(const cbmc::TrialMol& trialMo) const;
   
   //calculate correction term after swap move
   virtual double SwapCorrection(const cbmc::TrialMol& trialMol) const; 

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

   //update the hx,y,z hsqr and prefact
   virtual void UpdateRecipVec(uint box);

   //restore cosMol and sinMol
   virtual void RestoreMol(int molIndex);

   //update sinMol and cosMol
   virtual void exgMolCache();

};


#endif /*NOEWALD_H*/
