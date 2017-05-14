/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Ewald.h"
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

Ewald::Ewald(StaticVals const& stat, System & sys) : EwaldCached(stat, sys) {}

void Ewald::Init()
{
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

   electrostatic = forcefield.electrostatic;
   ewald = forcefield.ewald; 
   alpha = forcefield.alpha;
   recip_rcut = forcefield.recip_rcut;
   recip_rcut_Sq = recip_rcut * recip_rcut;
   SetNull();
   AllocMem();
   //initialize K vectors and reciprocate terms
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      RecipInit(b, currentAxes);
      BoxReciprocalSetup(b, currentCoords);
      SetRecipRef(b);
      printf("box: %d, RecipVectors: %d, kmax: %d\n", b, imageSize[b],
	     kmax[b]);
   }      
}


void Ewald::AllocMem()
{  
  //get size of image using defined Kmax  
  //Allocate Memory

  kmax = new uint[BOX_TOTAL];
  imageSize = new uint[BOX_TOTAL];
  imageSizeRef = new uint[BOX_TOTAL];
  sumRnew = new double*[BOX_TOTAL];
  sumInew = new double*[BOX_TOTAL];
  sumRref = new double*[BOX_TOTAL];
  sumIref = new double*[BOX_TOTAL];
  kx = new double*[BOX_TOTAL];
  ky = new double*[BOX_TOTAL];
  kz = new double*[BOX_TOTAL];
  hsqr = new double*[BOX_TOTAL];
  prefact = new double*[BOX_TOTAL];
  kxRef = new double*[BOX_TOTAL];
  kyRef = new double*[BOX_TOTAL];
  kzRef = new double*[BOX_TOTAL];
  hsqrRef = new double*[BOX_TOTAL];
  prefactRef = new double*[BOX_TOTAL];
     
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
     kx[b] = new double[imageTotal];
     ky[b] = new double[imageTotal];
     kz[b] = new double[imageTotal];
     hsqr[b] = new double[imageTotal];
     prefact[b] = new double[imageTotal];
     kxRef[b] = new double[imageTotal];
     kyRef[b] = new double[imageTotal];
     kzRef[b] = new double[imageTotal];
     hsqrRef[b] = new double[imageTotal];
     prefactRef[b] = new double[imageTotal];
     sumRnew[b] = new double[imageTotal];
     sumInew[b] = new double[imageTotal];
     sumRref[b] = new double[imageTotal];
     sumIref[b] = new double[imageTotal];     
  }      
}


void Ewald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
   uint counter = 0;
   int x, y, z, nky_max, nky_min, nkz_max, nkz_min;   
   double ksqr;
   double alpsqr4 = 1.0 / (4.0 * alpha * alpha);
   XYZ constValue = boxAxes.axis.Get(box);
   constValue.Inverse();
   constValue *= 2 * M_PI;
   double vol = boxAxes.volume[box] / (4 * M_PI);
   kmax[box] = int(recip_rcut * boxAxes.axis.BoxSize(box) / (2 * M_PI)) + 1;
 
   for (x = 0; x <= kmax[box]; x++)
   {
      nky_max = sqrt(pow(kmax[box], 2) - pow(x, 2));
      nky_min = -nky_max;
      if (x == 0.0)
      { 
	 nky_min = 0;
      }
      for (y = nky_min; y <= nky_max; y++)
      {
	 nkz_max = sqrt(pow(kmax[box], 2) - pow(x, 2) - pow(y, 2));
	 nkz_min = -nkz_max;
	 if (x == 0.0 && y == 0.0)
         { 
	    nkz_min = 1;
	 }
	 for (z = nkz_min; z <= nkz_max; z++)
         {
	   ksqr = pow((constValue.x * x), 2) + pow((constValue.y * y), 2) +
	     pow ((constValue.z * z), 2);
	    
	    if (ksqr < recip_rcut_Sq)
	    {
	       kx[box][counter] = constValue.x * x;
	       ky[box][counter] = constValue.y * y;
	       kz[box][counter] = constValue.z * z;
	       hsqr[box][counter] = ksqr;
	       prefact[box][counter] = num::qqFact * exp(-ksqr * alpsqr4)/
		 (ksqr * vol);
	       counter++;
	    }
	 }
      }
   }

   imageSize[box] = counter;
   
   if (counter > imageTotal)
   {
     std::cout<< "Error: Number of reciprocate vectors is greater than initialized vector size." << std::endl;  
     exit(EXIT_FAILURE);
   }
}


//calculate reciprocate term for a box
void Ewald::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
   uint i, j, m;
   double dotProduct = 0.0;
   double sumReal = 0.0;
   double sumImaginary = 0.0;

   if (box < BOXES_WITH_U_NB)
   {
      MoleculeLookup::box_iterator end = molLookup.BoxEnd(box);	 
      MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box);

#ifdef _OPENMP      
#pragma omp parallel default(shared) 
#endif
      {
	 std::memset(sumRnew[box], 0.0, sizeof(double) * imageSize[box]);
	 std::memset(sumInew[box], 0.0, sizeof(double) * imageSize[box]);
      }

      while (thisMol !=end)
      {	 
	 MoleculeKind const& thisKind = mols.GetKind(*thisMol);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, j, dotProduct, sumReal, sumImaginary)
#endif 
	 for (i = 0; i < imageSize[box]; i++)
	 {
	    sumReal = 0.0;
	    sumImaginary = 0.0;
 
	    for (j = 0; j < thisKind.NumAtoms(); j++)
	    {
	      dotProduct = currentAxes.DotProduct(mols.MolStart(*thisMol) + j,
						  kx[box][i], ky[box][i],
						  kz[box][i], molCoords, box);

	      sumReal += (thisKind.AtomCharge(j) * cos(dotProduct));
	      sumImaginary += (thisKind.AtomCharge(j) * sin(dotProduct));
	    }
	    sumRnew[box][i] += sumReal;
	    sumInew[box][i] += sumImaginary;
	 }
	 thisMol++;
      }
   }
}


//calculate reciprocate term for a box
double Ewald::BoxReciprocal(uint box) const
{
   uint i;
   double energyRecip = 0.0; 

   if (box < BOXES_WITH_U_NB)
   {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i) reduction(+:energyRecip)
#endif
      for (i = 0; i < imageSize[box]; i++)
      {
	energyRecip += (( sumRnew[box][i] * sumRnew[box][i] +
			  sumInew[box][i] * sumInew[box][i]) *
			prefact[box][i]);	
      }
   }

   return energyRecip; 
}


//calculate reciprocate term for displacement and rotation move
double Ewald::MolReciprocal(XYZArray const& molCoords,
			    const uint molIndex, const uint box,
			    XYZ const*const newCOM)
{
   double energyRecipNew = 0.0; 
   double energyRecipOld = 0.0;
   
   if (box < BOXES_WITH_U_NB)
   {
      MoleculeKind const& thisKind = mols.GetKind(molIndex);
      uint length = thisKind.NumAtoms();
      uint startAtom = mols.MolStart(molIndex);
      uint i, p, atom;
      double sumRealNew, sumImaginaryNew, dotProductNew, dotProductOld,
	sumRealOld, sumImaginaryOld;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, p, atom, sumRealNew, sumImaginaryNew, sumRealOld, sumImaginaryOld, dotProductNew, dotProductOld) reduction(+:energyRecipNew, energyRecipOld)
#endif
      for (i = 0; i < imageSizeRef[box]; i++)
      { 
	 sumRealNew = 0.0;
	 sumImaginaryNew = 0.0;
	 dotProductNew = 0.0;
	 dotProductOld = 0.0;
	 sumRealOld = 0.0;
	 sumImaginaryOld = 0.0;  
	    
	 for (p = 0; p < length; ++p)
	 {
	    atom = startAtom + p;
	    dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
						   kyRef[box][i], kzRef[box][i],
						   molCoords, box);

	    dotProductOld = currentAxes.DotProduct(atom, kxRef[box][i],
						   kyRef[box][i],kzRef[box][i],
						   currentCoords, box);
	    
	    sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
	    sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));

	    sumRealOld += (thisKind.AtomCharge(p) * cos(dotProductOld));
	    sumImaginaryOld += (thisKind.AtomCharge(p) * sin(dotProductOld));
	 }
	 
	 sumRnew[box][i] = sumRref[box][i] - sumRealOld + sumRealNew;
	 sumInew[box][i] = sumIref[box][i] - sumImaginaryOld + sumImaginaryNew;
	 
	 energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
			    * sumInew[box][i]) * prefactRef[box][i];	 
      }
   }

   return energyRecipNew - sysPotRef.boxEnergy[box].recip; 
}



//calculate reciprocate term in destination box for swap move
double Ewald::SwapDestRecip(const cbmc::TrialMol &newMol,
			    const uint box, const int sourceBox,
			    const int molIndex) 
{
   double energyRecipNew = 0.0; 
   double energyRecipOld = 0.0; 

    
   if (box < BOXES_WITH_U_NB)
   {
      uint p, i, length;
      MoleculeKind const& thisKind = newMol.GetKind();
      XYZArray molCoords = newMol.GetCoords();
      double dotProductNew, sumRealNew, sumImaginaryNew;
      length = thisKind.NumAtoms();

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, p, dotProductNew, sumRealNew, sumImaginaryNew) reduction(+:energyRecipNew) 
#endif
      for (i = 0; i < imageSizeRef[box]; i++)
      {
	 sumRealNew = 0.0;
	 sumImaginaryNew = 0.0;
	 dotProductNew = 0.0;  
	
	 for (p = 0; p < length; ++p)
	 {
	    dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
						   kyRef[box][i], kzRef[box][i],
						   molCoords, box);

	    sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
	    sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));
	 }

	 //sumRealNew;
	 sumRnew[box][i] = sumRref[box][i] + sumRealNew;
	 //sumImaginaryNew;
	 sumInew[box][i] = sumIref[box][i] + sumImaginaryNew;
   
	 energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
			    * sumInew[box][i]) * prefactRef[box][i];
      }

      energyRecipOld = sysPotRef.boxEnergy[box].recip;
   }

   return energyRecipNew - energyRecipOld;
}


//calculate reciprocate term in source box for swap move
double Ewald::SwapSourceRecip(const cbmc::TrialMol &oldMol,
			      const uint box, const int molIndex) 
{
   double energyRecipNew = 0.0; 
   double energyRecipOld = 0.0; 
   
   if (box < BOXES_WITH_U_NB)
   {
      uint i, p;
      double sumRealNew, sumImaginaryNew, dotProductNew;
      MoleculeKind const& thisKind = oldMol.GetKind();
      XYZArray molCoords = oldMol.GetCoords();
      uint length = thisKind.NumAtoms();

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, p, dotProductNew, sumRealNew, sumImaginaryNew) reduction(+:energyRecipNew)
#endif 
      for (i = 0; i < imageSizeRef[box]; i++)
      { 
	 sumRealNew = 0.0;
	 sumImaginaryNew = 0.0;
	 dotProductNew = 0.0;

	 for (p = 0; p < length; ++p)
	 {
	    dotProductNew = currentAxes.DotProduct(p, kxRef[box][i],
						   kyRef[box][i], kzRef[box][i],
						   molCoords, box);
	    
	    sumRealNew += (thisKind.AtomCharge(p) * cos(dotProductNew));
	    sumImaginaryNew += (thisKind.AtomCharge(p) * sin(dotProductNew));

	 }	
	 sumRnew[box][i] = sumRref[box][i] - sumRealNew;
	 sumInew[box][i] = sumIref[box][i] - sumImaginaryNew;
	
	 energyRecipNew += (sumRnew[box][i] * sumRnew[box][i] + sumInew[box][i]
			    * sumInew[box][i]) * prefactRef[box][i];	 
      }

      energyRecipOld = sysPotRef.boxEnergy[box].recip;
   }
   return energyRecipNew - energyRecipOld;
}


//restore cosMol and sinMol
void Ewald::RestoreMol(int molIndex)
{
   return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Ewald::exgMolCache()
{
   return;
}
