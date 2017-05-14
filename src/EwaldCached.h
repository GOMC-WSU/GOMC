/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EWALDCACHED_H
#define EWALDCACHED_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#include <omp.h>

#include "Molecules.h"
#include "Forcefield.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "MoleculeKind.h"
#include "TrialMol.h"


//
//    Calculating Electrostatic calculation with caching Fourier terms.
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

namespace cbmc { class TrialMol; }
namespace config_setup{ class SystemVals; }

class EwaldCached 
{
   public:

   EwaldCached(StaticVals const& stat, System & sys);
    
   ~EwaldCached();

   virtual void Init();

   virtual void AllocMem();
   void SetNull();

   //return size of image with defined Kmax value
   //uint GetImageSize();

   //initiliazie term used for ewald calculation
   virtual void RecipInit(uint box, BoxDimensions const& boxAxes);
   virtual void RecipCountInit(uint box, BoxDimensions const& boxAxes);
   
   //calculate self term for a box
   virtual double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

   //setup reciprocate term for a box
   virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

   //calculate reciprocate energy term for a box
   virtual double BoxReciprocal(uint box) const;

   //calculate reciprocate force term for a box
   virtual Virial ForceReciprocal(Virial& virial, uint box) const;


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

   //update the hx,y,z hsqr and prefact
   virtual void UpdateRecipVec(uint box);

   //restore cosMol and sinMol
   virtual void RestoreMol(int molIndex);

   //update sinMol and cosMol
   virtual void exgMolCache();

   uint findLargeImage();

   protected: 
   
   const Forcefield& forcefield;
   const Molecules& mols;
   const Coordinates& currentCoords;
   const MoleculeLookup& molLookup;
   const BoxDimensions& currentAxes;
   const COM& currentCOM;
   const SystemPotential &sysPotRef;
   

   bool electrostatic, ewald;
   double alpha; 
   double recip_rcut, recip_rcut_Sq;
   uint *imageSize;
   uint *imageSizeRef;
   //const uint imageTotal = GetImageSize();
   const int imageTotal;
   uint memoryAllocation;
   uint imageLarge;
   uint *kmax;
   double **sumRnew; //cosine serries
   double **sumInew; //sine serries
   double **sumRref;
   double **sumIref;
   double *cosMolRestore; //cos()*charge
   double *sinMolRestore; //sin()*charge
   double **cosMolRef;
   double **sinMolRef;
   double **cosMolBoxRecip;
   double **sinMolBoxRecip;
   double **kx, **kxRef;
   double **ky, **kyRef; 
   double **kz, **kzRef;
   double **hsqr, **hsqrRef;
   double **prefact, **prefactRef;
   

   std::vector<int> particleKind;
   std::vector<int> particleMol;
   std::vector<double> particleCharge;
   
};


#endif /*EWALDCACHED_H*/
