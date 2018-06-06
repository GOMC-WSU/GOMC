#ifndef MOLECULEEXCHANGE1_H
#define MOLECULEEXCHANGE1_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include <cmath>

using std::vector;
using namespace geom;

// MEMC-1 Move:
//
// Swapping one Large molecule with one or more small molecules in dense phase 
// and vice versa.
// Sub-Volume location and orientation is random.

class MoleculeExchange1 : public MoveBase
{
 public:

   MoleculeExchange1(System &sys, StaticVals const& statV) :
   ffRef(statV.forcefield), molLookRef(sys.molLookupRef), MoveBase(sys, statV),
   cavity(statV.memcVal.subVol), cavA(3), invCavA(3), kindS(-1), kindL(-1),
   perAdjust(statV.GetPerAdjust())
   {
     enableID = statV.memcVal.enable;
     largeBB[0] = -1;
     largeBB[1] = -1;

     if(enableID) {
       if(molLookRef.GetNumCanSwapKind() < 2) {
	 std::cout << "Error: MEMC move cannot be applied to pure systems or"<<
	   " systems, where only one molecule type is allowed to be swapped.\n";
	 exit(EXIT_FAILURE);
       }

       if(cavity.x >= cavity.y)
         cavity.y = cavity.x;
       else
         cavity.x = cavity.y;
         
       volCav = cavity.x * cavity.y * cavity.z;
       exchangeRatio = statV.memcVal.exchangeRatio;
         
       for(uint k = 0; k < molLookRef.GetNumCanSwapKind(); k++) {
         if(molRef.kinds[k].name == statV.memcVal.largeKind) {
           kindL = molLookRef.GetCanSwapKind(k);
         } else if(molRef.kinds[k].name == statV.memcVal.smallKind) {
	   kindS = molLookRef.GetCanSwapKind(k);
         }
       }

       if(kindS == -1) {
	 printf("Error: Residue name %s was not found in PDB file as small molecule kind to be exchanged.\n", statV.memcVal.smallKind.c_str());
	 exit(EXIT_FAILURE);
       }

       if(kindL == -1) {
	 printf("Error: Residue name %s was not found in PDB file as large molecule kind to be exchanged.\n", statV.memcVal.largeKind.c_str());
	 exit(EXIT_FAILURE);
       }
         
       for(uint i = 0; i < molRef.kinds[kindL].NumAtoms(); i++) {
	 if(molRef.kinds[kindL].atomNames[i] == statV.memcVal.largeBBAtom1) {
	   largeBB[0] = i;
	 }
	 else if(molRef.kinds[kindL].atomNames[i] ==statV.memcVal.largeBBAtom2){
	   largeBB[1] = i;
	 }
       }
         
       for(uint i = 0; i < 2; i++) {
	 if(largeBB[i] == -1) {
	   printf("Error: Atom name %s or %s was not found in %s residue.\n",
		  statV.memcVal.largeBBAtom1.c_str(),
		  statV.memcVal.largeBBAtom2.c_str(),
		  statV.memcVal.largeKind.c_str());
	   exit(EXIT_FAILURE);
	 }
       }
     }

     //checking the acceptance statistic for each kind
     counter = 0;
     molInCavCount = 0;
     lastAccept = 0.0;
     exDiff = 1;
   }

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);

 private:

   void AdjustExRatio();
   void ShiftMol(const bool A, const uint n, const uint from, const uint to);
   void RecoverMol(const bool A, const uint n, const uint from, const uint to);
   uint PickMolInCav();
   uint ReplaceMolecule();
   void CalcTc();
   double GetCoeff() const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);

   bool insertL, enableID;
   uint largeBB[2];
   uint sourceBox, destBox;
   uint perAdjust, molInCavCount, counter;
   uint numInCavA, numInCavB, kindS, kindL, totMolInCav;
   vector<uint> pStartA, pLenA, pStartB, pLenB;
   vector<uint> molIndexA, kindIndexA, molIndexB, kindIndexB;
   vector< vector<uint> > molInCav;
   vector<cbmc::TrialMol> oldMolA, newMolA, oldMolB, newMolB;

   int exDiff, exchangeRatio;
   double volCav, lastAccept;
   double numTypeASource, numTypeBSource, numTypeADest, numTypeBDest;
   XYZ center, cavity;
   XYZArray cavA, invCavA;
   double W_tc, W_recip;
   double correct_oldA, correct_newA, self_oldA, self_newA;
   double correct_oldB, correct_newB, self_oldB, self_newB;
   double recipDest, recipSource;
   Intermolecular tcNew[BOX_TOTAL];
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline void MoleculeExchange1::AdjustExRatio()
{
  if(((counter + 1) % perAdjust) == 0)
  {
    uint exMax = ceil((float)molInCavCount / (float)perAdjust);
    uint exMin = ceil((float)exMax / 2.0);
    if(exMin == 0)
      exMin = 1;

    subPick = mv::GetMoveSubIndex(mv::MEMC, sourceBox);
    double currAccept = moveSetRef.GetAccept(subPick);
    if(abs(currAccept - lastAccept) > 0.05 * currAccept)
    {
      if(currAccept > lastAccept)
      {
	exchangeRatio += exDiff;
      }
      else
      {
	exDiff *= -1;
	exchangeRatio += exDiff;
      }
      lastAccept = currAccept;
      if(exchangeRatio < exMin)
	exchangeRatio  = exMin;
      if(exchangeRatio > exMax)
	exchangeRatio = exMax;
    }
    molInCavCount = 0;
    counter = 0;
    printf("Average Mol In Cavity: %d. Exchange Ratio: %d \n", exMax,
	   exchangeRatio);
  }
}


inline uint MoleculeExchange1::PickMolInCav()
{
   uint state = mv::fail_state::NO_FAIL;
   //pick a random location in dense phase
   XYZ axis = boxDimRef.GetAxis(sourceBox);
   XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
   //Use to shift the new insterted molecule
   center = temp;
   //Pick random vector anad find two vectors that are perpendicular to V1
   SetBasis(cavA, prng.RandomUnitVect());
   //Calculate inverse matrix for cav here Inv = transpose
   TransposeMatrix(invCavA, cavA);

   //Find the molecule kind 0 in the cavity
   if(calcEnRef.FindMolInCavity(molInCav, center, cavity, invCavA,
				sourceBox, kindS, exchangeRatio))
   {
     molIndexA.clear();
     kindIndexA.clear();
     molIndexB.clear();
     kindIndexB.clear();
     //Find the exchangeRatio number of molecules kind 0 in cavity
     numInCavA = exchangeRatio;
     totMolInCav = molInCav[kindS].size();
     for(uint n = 0; n < numInCavA; n++)
     {
       //pick random exchangeRatio number of kindS in cavity
       uint picked = prng.randIntExc(molInCav[kindS].size());
       molIndexA.push_back(molInCav[kindS][picked]);
       kindIndexA.push_back(molRef.GetMolKind(molIndexA[n]));
       molInCav[kindS].erase(molInCav[kindS].begin() + picked);
     } 

     //pick a molecule from Large kind in destBox
     numInCavB = 1;
     state = prng.PickMol(kindL, kindIndexB, molIndexB, numInCavB, destBox);
   }
   else
   {
     //reject the move
     state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
   }

   molInCavCount += totMolInCav;
   counter++;

   return state;
}


inline uint MoleculeExchange1::ReplaceMolecule()
{
   uint state = mv::fail_state::NO_FAIL;
   molIndexA.clear();
   kindIndexA.clear();
   molIndexB.clear();
   kindIndexB.clear();
   numInCavA = 1;
   numInCavB = exchangeRatio;
   //pick a random molecule of Large kind in dens box
   state = prng.PickMol(kindL, kindIndexA, molIndexA, numInCavA, sourceBox);

   if(state == mv::fail_state::NO_FAIL)
   {
     //Set the V1 to the backbone of the large molecule
     if(molRef.NumAtoms(kindL) == 1)
     {
       SetBasis(cavA, prng.RandomUnitVect());
     }
     else
     {
       SetBasis(cavA, boxDimRef.MinImage(coordCurrRef.Difference(largeBB[0], largeBB[1]), sourceBox));
     }
     //Calculate inverse matrix for cav. Here Inv = Transpose 
     TransposeMatrix(invCavA, cavA);
     //Use to shift to the COM of new molecule
     center = comCurrRef.Get(molIndexA[0]);
     //find how many of KindS exist in this center
     calcEnRef.FindMolInCavity(molInCav, center, cavity, invCavA, sourceBox,
			       kindS, exchangeRatio);
     totMolInCav = molInCav[kindS].size();
     //pick exchangeRatio number of Small molecule from dest box
     state = prng.PickMol(kindS, kindIndexB, molIndexB, numInCavB, destBox);
   }
   return state;
}

inline uint MoleculeExchange1::GetBoxPairAndMol(const double subDraw,
						const double movPerc)
{
   uint state = mv::fail_state::NO_FAIL; 
   //deside to insert or remove the big molecule
   prng.PickBool(insertL, subDraw, movPerc);
   
#if ENSEMBLE == GEMC
   double density;
   double maxDens = 0.0;
   uint densB;
   //choose the sourceBox to be the dense phase
   for(uint b = 0; b < BOX_TOTAL; b++)
   {
     density = 0.0;
     for(uint k = 0; k < molLookRef.GetNumKind(); k++)
     {
       density += molLookRef.NumKindInBox(k, b) * boxDimRef.volInv[b] *
	 molRef.kinds[k].molMass;
     }
     if(density > maxDens)
     {
       maxDens = density;
       densB = b;
     }
   }

   //Pick box in dense phase
   sourceBox = densB; 
   //Pick the destination box
   prng.SetOtherBox(destBox, sourceBox);
   //prng.PickBoxPair(sourceBox, destBox, subDraw, movPerc);

#elif ENSEMBLE == GCMC
   sourceBox = 0;
   destBox = 1;
#endif

   //adjust exchange rate based on number of small kind in cavity
   //AdjustExRatio();

   if(insertL)
   {
     state = PickMolInCav();
   }
   else
   {
     state = ReplaceMolecule();
   }  
   
   if(state == mv::fail_state::NO_FAIL)
   {
     pStartA.clear();
     pStartB.clear();
     pStartA.resize(numInCavA);
     pStartB.resize(numInCavB);
     pLenA.clear();
     pLenB.clear();
     pLenA.resize(numInCavA);
     pLenB.resize(numInCavB);

     for(uint n = 0; n < numInCavA; n++)
     {
       pStartA[n] = pLenA[n] = 0;
       molRef.GetRangeStartLength(pStartA[n], pLenA[n], molIndexA[n]);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       pStartB[n] = pLenB[n] = 0;
       molRef.GetRangeStartLength(pStartB[n], pLenB[n], molIndexB[n]);
     }
   }

   return state;
}


inline uint MoleculeExchange1::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   if(state == mv::fail_state::NO_FAIL)
   {
     newMolA.clear();
     oldMolA.clear();
     newMolB.clear();
     oldMolB.clear();

     numTypeASource =(double)(molLookRef.NumKindInBox(kindIndexA[0],sourceBox));
     numTypeADest = (double)(molLookRef.NumKindInBox(kindIndexA[0], destBox));
     numTypeBSource =(double)(molLookRef.NumKindInBox(kindIndexB[0],sourceBox));
     numTypeBDest = (double)(molLookRef.NumKindInBox(kindIndexB[0], destBox));

     //transfering type A from source to dest
     for(uint n = 0; n < numInCavA; n++)
     {
       newMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					destBox));
       oldMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
					sourceBox));
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       //transfering type B from dest to source
       newMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					sourceBox));
       oldMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
					destBox));
     }

     //set the old coordinate after unwrap them
     for(uint n = 0; n < numInCavA; n++)
     {
       XYZArray molA(pLenA[n]);
       coordCurrRef.CopyRange(molA, pStartA[n], 0, pLenA[n]);
       boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA[n]));
       oldMolA[n].SetCoords(molA, 0);
       //set coordinate of moleA to newMolA, later it will shift to center
       newMolA[n].SetCoords(molA, 0); 
       //copy cavA matrix to slant the old trial of molA
       oldMolA[n].SetCavMatrix(cavA);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       XYZArray molB(pLenB[n]);     
       coordCurrRef.CopyRange(molB, pStartB[n], 0, pLenB[n]);
       boxDimRef.UnwrapPBC(molB, destBox, comCurrRef.Get(molIndexB[n]));
       oldMolB[n].SetCoords(molB, 0);
       //set coordinate of moleB to newMolB, later it will shift
       newMolB[n].SetCoords(molB, 0);
       //copy cavA matrix to slant the new trial of molB
       newMolB[n].SetCavMatrix(cavA);
     }

     for(uint n = 0; n < numInCavB; n++)
     {
       //SetSeed(has cavity, COM is fixed, rotate around Backbone)
       if(insertL)
       {
	 //Inserting Lmol from destBox to the center of cavity in sourceBox
	 newMolB[n].SetSeed(center, cavity, true, true, true);
	 //perform rotational trial move in destBox for L oldMol
	 oldMolB[n].SetSeed(false, false, false);
       }
       else
       {
	 //Inserting S mol from destBox to the cavity in sourceBox
	 newMolB[n].SetSeed(center, cavity, true, false, false);
	 //perform trial move in destBox for S oldMol
	 oldMolB[n].SetSeed(false, false, false);
       }
     }

     for(uint n = 0; n < numInCavA; n++)
     {
       if(insertL)
       {
	 //Inserting S mol from sourceBox to destBox
	 newMolA[n].SetSeed(false, false, false);
	 ////perform trial move in cavity in sourceBox for S oldMol
	 oldMolA[n].SetSeed(center, cavity, true, false, false);
       }
       else
       {
	 //Inserting L mol from sourceBox to destBox
	 newMolA[n].SetSeed(false, false, false);
	 //perform rotational trial move on COM for L oldMol
	 oldMolA[n].SetSeed(center, cavity, true, true, true);
       }
     }
   }

   return state;
}


inline uint MoleculeExchange1::Transform()
{
  CalcTc();

  //Calc Old energy and delete A from source
  for(uint n = 0; n < numInCavA; n++)
  {
    cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexA[n]].BuildIDOld(oldMolA[n], molIndexA[n]);
    //Add bonded energy because we dont considered in DCRotate.cpp
    calcEnRef.MoleculeIntra(oldMolA[n], molIndexA[n]);
  }
  
  //Calc old energy and delete B from destBox
  for(uint n = 0; n < numInCavB; n++)
  {
    cellList.RemoveMol(molIndexB[n], destBox, coordCurrRef);
    molRef.kinds[kindIndexB[n]].BuildIDOld(oldMolB[n], molIndexB[n]);
    //Add bonded energy because we dont considered in DCRotate.cpp
    calcEnRef.MoleculeIntra(oldMolB[n], molIndexB[n]);
  }
  
  //Insert A to destBox
  for(uint n = 0; n < numInCavA; n++)
  {
    molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
    ShiftMol(true, n, sourceBox, destBox);
    cellList.AddMol(molIndexA[n], destBox, coordCurrRef);
    //Add bonded energy because we dont considered in DCRotate.cpp
    calcEnRef.MoleculeIntra(newMolA[n], molIndexA[n]);
  }

  //Insert B in sourceBox
  for(uint n = 0; n < numInCavB; n++)
  {
    molRef.kinds[kindIndexB[n]].BuildIDNew(newMolB[n], molIndexB[n]);
    ShiftMol(false, n, destBox, sourceBox);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);      
    //Add bonded energy because we dont considered in DCRotate.cpp
    calcEnRef.MoleculeIntra(newMolB[n], molIndexB[n]); 
  }
  
  return mv::fail_state::NO_FAIL;
}

inline void MoleculeExchange1::CalcTc()
{
  W_tc = 1.0;
  if (ffRef.useLRC)
  {
    double delTC = 0.0;
    for (uint b = 0; b < BOX_TOTAL; ++b)
    {
      uint kCount[molRef.kindsCount];
      for (uint k = 0; k < molRef.kindsCount; ++k)
      {
	kCount[k] = molLookRef.NumKindInBox(k, b);
      }

      if (b == sourceBox)
      {
	kCount[kindIndexA[0]] -= numInCavA;
	kCount[kindIndexB[0]] += numInCavB;	   
      }
      else if (b == destBox)
      {
	kCount[kindIndexA[0]] += numInCavA;
	kCount[kindIndexB[0]] -= numInCavB;
      }
      tcNew[b].energy = calcEnRef.EnergyCorrection(b, kCount);
      delTC += tcNew[b].energy - sysPotRef.boxEnergy[b].tc;
    }
    W_tc = exp(-1.0 * ffRef.beta * delTC); 
  }
}
inline void MoleculeExchange1::CalcEn()
{   
   W_recip = 1.0;
   correct_oldA = 0.0, correct_newA = 0.0;
   self_oldA = 0.0, self_newA = 0.0;
   correct_oldB = 0.0, correct_newB = 0.0;
   self_oldB = 0.0, self_newB = 0.0;
   recipDest = 0.0, recipSource = 0.0;

   for(uint n = 0; n < numInCavA; n++)
   {
      correct_newA += calcEwald->SwapCorrection(newMolA[n]);
      correct_oldA += calcEwald->SwapCorrection(oldMolA[n]);
      self_newA += calcEwald->SwapSelf(newMolA[n]);
      self_oldA += calcEwald->SwapSelf(oldMolA[n]);
   }
   recipDest = calcEwald->SwapRecip(newMolA, oldMolB);

   for(uint n = 0; n < numInCavB; n++)
   {
     correct_newB += calcEwald->SwapCorrection(newMolB[n]);
     correct_oldB += calcEwald->SwapCorrection(oldMolB[n]);
     self_newB += calcEwald->SwapSelf(newMolB[n]);
     self_oldB += calcEwald->SwapSelf(oldMolB[n]);
   }
   recipSource = calcEwald->SwapRecip(newMolB, oldMolA);

   //need to contribute the self and correction energy 
   W_recip = exp(-1.0 * ffRef.beta * (recipSource + recipDest +
				      correct_newA - correct_oldA +
				      correct_newB - correct_oldB +
				      self_newA - self_oldA +
				      self_newB - self_oldB));
   
}

inline double MoleculeExchange1::GetCoeff() const
{
  double volSource = boxDimRef.volume[sourceBox];
  double volDest = boxDimRef.volume[destBox];
#if ENSEMBLE == GEMC
  if(insertL)
  {
    //kindA is the small molecule
    double ratioF =  num::Factorial(totMolInCav) /
      (num::Factorial(totMolInCav - exchangeRatio) *
       num:: Factorial(numTypeADest, exchangeRatio));

    double ratioV = (volSource / volDest) * pow(volDest / volCav, exchangeRatio);
    
    return ratioF * ratioV  * numTypeBDest / (numTypeBSource + 1.0);
  }
  else
  {
    //kindA is the big molecule
    double ratioF = num::Factorial(totMolInCav) *
      num::Factorial(numTypeBDest - exchangeRatio, exchangeRatio ) /
      num::Factorial(totMolInCav + exchangeRatio);

    double ratioV = (volDest / volSource) * pow(volCav / volDest, exchangeRatio);

    return ratioF * ratioV  * numTypeASource / (numTypeADest + 1.0);
  }
#elif ENSEMBLE == GCMC
  if(ffRef.isFugacity)
  {
    double delA = molRef.kinds[kindIndexA[0]].chemPot * numInCavA;
    double insB = molRef.kinds[kindIndexB[0]].chemPot * numInCavB;
    if(insertL)
    {
      //Insert Large molecule
      double ratioF = num::Factorial(totMolInCav) /
	num::Factorial(totMolInCav - exchangeRatio);

      double ratioV = volSource / pow(volCav, exchangeRatio);
      return (insB / delA) * ratioF * ratioV / (numTypeBSource + 1.0);
    }
    else
    {
      //Delete Large Molecule
      double ratioF = num::Factorial(totMolInCav) /
	num::Factorial(totMolInCav + exchangeRatio);

      double ratioV = pow(volCav, exchangeRatio) / volSource;
      return (insB / delA) * ratioF * ratioV * numTypeASource;
    }
  }
  else
  {
    double delA = (-BETA * molRef.kinds[kindIndexA[0]].chemPot * numInCavA);
    double insB = (BETA * molRef.kinds[kindIndexB[0]].chemPot * numInCavB);
    if(insertL)
    {
      // Insert Large molecule
      double ratioF = num::Factorial(totMolInCav) /
	num::Factorial(totMolInCav - exchangeRatio);

      double ratioV = volSource / pow(volCav, exchangeRatio);
      return exp(delA + insB) * ratioF * ratioV / (numTypeBSource + 1.0);
    }
    else
    {
      //Delete Large molecule
      double ratioF = num::Factorial(totMolInCav) /
	num::Factorial(totMolInCav + exchangeRatio);

      double ratioV = pow(volCav, exchangeRatio) / volSource;
      return exp(delA + insB) * ratioF * ratioV * numTypeASource;
    }
  }
#endif
}

inline void MoleculeExchange1::ShiftMol(const bool A, const uint n,
				       const uint from, const uint to)
{
  if(A)
  {
    //Add type A to dest box
    newMolA[n].GetCoords().CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], to);
    molLookRef.ShiftMolBox(molIndexA[n], from, to, kindIndexA[n]);
  }
  else
  {
    //Add type B to source box
    newMolB[n].GetCoords().CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], to);
    molLookRef.ShiftMolBox(molIndexB[n], from, to, kindIndexB[n]);
  }
}

inline void MoleculeExchange1::RecoverMol(const bool A, const uint n,
					 const uint from, const uint to)
{
  if(A)
  {
    XYZArray molA(pLenA[n]);
    oldMolA[n].GetCoords().CopyRange(molA, 0, 0, pLenA[n]);
    boxDimRef.WrapPBC(molA, to);

    molA.CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], to);
    molLookRef.ShiftMolBox(molIndexA[n], from, to, kindIndexA[n]);
  }
  else
  {
    XYZArray molB(pLenB[n]);
    oldMolB[n].GetCoords().CopyRange(molB, 0, 0, pLenB[n]);
    boxDimRef.WrapPBC(molB, to);

    molB.CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], to);
    molLookRef.ShiftMolBox(molIndexB[n], from, to, kindIndexB[n]);
  }
}

inline void MoleculeExchange1::Accept(const uint rejectState, const uint step)
{
   bool result;

   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();  
      double Wrat = W_tc * W_recip;

      for(uint n = 0; n < numInCavA; n++)
      {
	Wrat *= newMolA[n].GetWeight() / oldMolA[n].GetWeight();
      }

      for(uint n = 0; n < numInCavB; n++)
      {
	Wrat *= newMolB[n].GetWeight() / oldMolB[n].GetWeight();
      }

      result = prng() < molTransCoeff * Wrat;

      if(result)
      {
	 //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc = tcNew[sourceBox].energy;
         sysPotRef.boxEnergy[destBox].tc = tcNew[destBox].energy;

         //Add rest of energy.
	 for(uint n = 0; n < numInCavB; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] += newMolB[n].GetEnergy();
	   sysPotRef.boxEnergy[destBox] -= oldMolB[n].GetEnergy();
	 }

	 for(uint n = 0; n < numInCavA; n++)
	 {
	   sysPotRef.boxEnergy[sourceBox] -= oldMolA[n].GetEnergy();
	   sysPotRef.boxEnergy[destBox] += newMolA[n].GetEnergy();
	 }
	 

	 //Add Reciprocal energy
	 sysPotRef.boxEnergy[sourceBox].recip += recipSource;
	 sysPotRef.boxEnergy[destBox].recip += recipDest;	 
	 //Add correction energy
	 sysPotRef.boxEnergy[sourceBox].correction -= correct_oldA;
	 sysPotRef.boxEnergy[sourceBox].correction += correct_newB;
	 sysPotRef.boxEnergy[destBox].correction += correct_newA;
	 sysPotRef.boxEnergy[destBox].correction -= correct_oldB;	 
	 //Add self energy
	 sysPotRef.boxEnergy[sourceBox].self -= self_oldA;
	 sysPotRef.boxEnergy[sourceBox].self += self_newB;
	 sysPotRef.boxEnergy[destBox].self += self_newA;
	 sysPotRef.boxEnergy[destBox].self -= self_oldB;
	 
	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    calcEwald->UpdateRecip(b);
	 }

	 //molA and molB already transfered to destBox and added to cellist

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
	//transfer molA from destBox to source
	for(uint n = 0; n < numInCavA; n++)
	{
	  cellList.RemoveMol(molIndexA[n], destBox, coordCurrRef);
	  RecoverMol(true, n, destBox, sourceBox);
	  cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
	}
	//transfer molB from sourceBox to dest
	for(uint n = 0; n < numInCavB; n++)
	{
	  cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
	  RecoverMol(false, n, sourceBox, destBox);
	  cellList.AddMol(molIndexB[n], destBox, coordCurrRef);
	}
      }
   }
   else  //else we didn't even try because we knew it would fail
      result = false;

#if ENSEMBLE == GEMC
   subPick = mv::GetMoveSubIndex(mv::MEMC, sourceBox);
   moveSetRef.Update(result, subPick, step);
   subPick = mv::GetMoveSubIndex(mv::MEMC, destBox);
   moveSetRef.Update(result, subPick, step);
#elif ENSEMBLE == GCMC
   subPick = mv::GetMoveSubIndex(mv::MEMC);
   moveSetRef.Update(result, subPick, step);
#endif
}

#endif

#endif        
