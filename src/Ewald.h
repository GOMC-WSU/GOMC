#ifndef EWALD_H
#define EWALD_H

#include "../lib/BasicTypes.h"
#include "EnergyTypes.h"
#include <vector>
#include <stdio.h>
#include <cstring>
//
//    Ewald.h
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocate part of ewald    
//
//    Developed by Y. Li and modified by Mohammad S. Barhaghi
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


class Ewald 
{
  //friend class CalculateEnergy;
   public:
    Ewald(StaticVals const& stat, System & sys);
    
    ~Ewald();

   void Init();

   void InitEwald();

   //return size of image with defined Kmax value
   //uint GetImageSize();

   //initiliazie term used for ewald calculation
   void RecipInit(uint box, BoxDimensions const& boxAxes);
   void RecipCountInit(uint box, BoxDimensions const& boxAxes);
   
   //calculate self term for a box
   double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

   //calculate reciprocate term for a box
   double BoxReciprocal(int box, XYZArray const& molCoords);

   //calculate correction term for a molecule
   double MolCorrection(uint molIndex, BoxDimensions const& boxAxes,
			uint box)const;

   //calculate reciprocate term for displacement and rotation move
   double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
			const uint box, XYZ const*const newCOM = NULL) ;	

   //calculate self term for CBMC algorithm
   void SwapSelf(double *self, uint molIndex, uint partIndex, int box, 
		 uint trials) const;
   
   //calculate correction term for linear molecule CBMC algorithm
   void SwapCorrection(double* energy, const cbmc::TrialMol& trialMol, 
		       XYZArray const& trialPos, const uint partIndex, 
		       const uint box, const uint trials) const; 

   //calculate correction term for branched molecule CBMC algorithm
   void SwapCorrection(double* energy, const cbmc::TrialMol& trialMol,
		       XYZArray *trialPos, const int pickedAtom, 
		       uint *partIndexArray, const uint box, const uint trials,
		       const uint PrevIndex, bool Prev) const;  

   //calculate correction term for old configuration
   double CorrectionOldMol(const cbmc::TrialMol& oldMol, const double distSq,
			   const uint i, const uint j) const;

   //calculate reciprocate term in destination box for swap move
   double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box, 
			const int sourceBox, const int molIndex);	

   //calculate reciprocate term in source box for swap move
   double SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box, const int molIndex);

   //back up reciprocate values
   void BackUpRecip( uint box)
   {  
     std::memcpy(sumRnew[box], sumRref[box], sizeof(double) * imageSize[box]);
     std::memcpy(sumInew[box], sumIref[box], sizeof(double) * imageSize[box]);
   }

   //update reciprocate values
   void UpdateRecip(uint box)
   {
     double *tempR, *tempI;
     tempR = sumRref[box];
     tempI = sumIref[box];
     sumRref[box] = sumRnew[box];
     sumIref[box] = sumInew[box];
     sumRnew[box] = tempR;
     sumInew[box] = tempI;
   }

   //restore cosMol and sinMol
   void RestoreMol(int molIndex)
   {
     double *tempCos, *tempSin;
     tempCos = cosMolRef[molIndex];
     tempSin = sinMolRef[molIndex];
     cosMolRef[molIndex] = cosMolRestore;
     sinMolRef[molIndex] = sinMolRestore;
     cosMolRestore = tempCos;
     sinMolRestore = tempSin;
    
     //    std::memcpy(cosMolRef[molIndex], cosMolRestore, sizeof(double)*imageLarge);
     //     std::memcpy(sinMolRef[molIndex], sinMolRestore, sizeof(double)*imageLarge);
   }
   
   int findLargeImage();
   void exgMolCache();

   private: 
   
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
   int *imageSize;
   //const uint imageTotal = GetImageSize();
   const int imageTotal;
   uint memoryAllocation;
   int imageLarge;
   int *kmax;
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
   double **kx;
   double **ky; 
   double **kz;
   double **prefact;
   

   std::vector<int> particleKind;
   std::vector<int> particleMol;
   std::vector<double> particleCharge;
   
};

#endif /*EWALD_H*/
