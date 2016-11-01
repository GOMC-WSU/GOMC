#include "NoEwald.h"
#include "EwaldCached.h"
#include "CalculateEnergy.h"        
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
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#include <omp.h>

//
//   
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald    
//
//    Developed by Y. Li and Mohammad S. Barhaghi
// 
//

using namespace geom;

NoEwald::NoEwald(StaticVals const& stat, System & sys) :
  EwaldCached(stat, sys) {}

void NoEwald::Init()
{

   electrostatic = forcefield.electrostatic;
   ewald = forcefield.ewald; 
   SetNull();
      
}


void NoEwald::AllocMem()
{  
   return;     
}


void NoEwald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
   return;   
}


//calculate reciprocate term for a box
void NoEwald::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
   return;   
}


//calculate reciprocate term for a box
double NoEwald::BoxReciprocal(uint box) const
{
   return 0.0; 
}


//calculate reciprocate term for displacement and rotation move
double NoEwald::MolReciprocal(XYZArray const& molCoords,
			      const uint molIndex, const uint box,
			      XYZ const*const newCOM)
{
   return 0.0; 
}


//calculate self term for a box
double NoEwald::BoxSelf(BoxDimensions const& boxAxes, uint box) const
{
   return 0.0;
}


//calculate correction term for a molecule
double NoEwald::MolCorrection(uint molIndex, BoxDimensions const& boxAxes,
			      uint box) const
{
   return 0.0;
}

//calculate reciprocate term in destination box for swap move
double NoEwald::SwapDestRecip(const cbmc::TrialMol &newMol,
			      const uint box, const int sourceBox,
			      const int molIndex) 
{
   return 0.0;
}


//calculate reciprocate term in source box for swap move
double NoEwald::SwapSourceRecip(const cbmc::TrialMol &oldMol,
				       const uint box, const int molIndex) 
{
   return 0.0;
}


//calculate self term for CBMC algorithm
void NoEwald::SwapSelf(double *self, uint molIndex, uint partIndex,
		       int box, uint trials) const
{
   return;
}

//calculate correction term for linear molecule CBMC algorithm
void NoEwald::SwapCorrection(double* energy,
			     const cbmc::TrialMol& trialMol, 
			     XYZArray const& trialPos,
			     const uint partIndex, const uint box,
			     const uint trials) const
{
   return;
}


//calculate correction term for branched molecule CBMC algorithm
void NoEwald::SwapCorrection(double* energy,
			     const cbmc::TrialMol& trialMol,
			     XYZArray *trialPos, const int pickedAtom, 
			     uint *partIndexArray, const uint box,
			     const uint trials,
			     const uint prevIndex, bool prev) const
{
   return;
}

double NoEwald::CorrectionOldMol(const cbmc::TrialMol& oldMol,
				 const double distSq, const uint i,
				 const uint j) const
{
   return 0.0;
}


//back up reciptocate value to Ref (will be called during initialization)
void NoEwald::SetRecipRef(uint box)
{
   return;
}

   //update reciprocate values
void NoEwald::UpdateRecip(uint box)
{
   return;
}

//restore cosMol and sinMol
void NoEwald::RestoreMol(int molIndex)
{
   return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::exgMolCache()
{
   return;
}


