/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOLCULETRANSFER_H
#define MOLCULETRANSFER_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "cbmc/TrialMol.h"

//#define DEBUG_MOVES

class MoleculeTransfer : public MoveBase
{
 public:

   MoleculeTransfer(System &sys, StaticVals const& statV) : 
      ffRef(statV.forcefield), molLookRef(sys.molLookupRef), 
      MoveBase(sys, statV) {}

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);
 private:
   double GetCoeff() const;
   uint GetBoxPairAndMol(const double subDraw, const double movPerc);
   MolPick molPick;
   uint sourceBox, destBox;
   uint pStart, pLen;
   uint molIndex, kindIndex;

   double W_tc, oldVirial;
   cbmc::TrialMol oldMol, newMol;
   Intermolecular tcLose, tcGain;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint MoleculeTransfer::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   uint state = prng.PickMolAndBoxPair(molIndex, kindIndex, sourceBox, destBox,
					  subDraw, movPerc);
 
   if ( state != mv::fail_state::NO_MOL_OF_KIND_IN_BOX)
   {
      pStart = pLen = 0;
      molRef.GetRangeStartLength(pStart, pLen, molIndex);
   }
   return state;
}

inline uint MoleculeTransfer::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxPairAndMol(subDraw, movPerc);
   newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
   oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
 // printf("kind index=%d, plen=%d\n", kindIndex,pLen);


   oldMol.SetCoords(coordCurrRef, pStart);
   W_tc = 1.0;
   return state;
}


inline uint MoleculeTransfer::Transform()
{
   oldVirial = calcEnRef.MoleculeVirial(molIndex, sourceBox);
   subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, sourceBox);
   molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);

  /* for (int i=0;i< 5;i++)
   {
   
   printf("x=%f,y=%f,z=%f\n", oldMol.tCoords.x[i],oldMol.tCoords.y[i],oldMol.tCoords.z[i] );
   }*/


 
   //printf("new energy=%f\n",newMol.en.inter);
   



   return mv::fail_state::NO_FAIL;
}

inline void MoleculeTransfer::CalcEn()
{
   if (ffRef.useLRC)
   {
       tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);

	 // printf("energy lost=%f\n",tcLose.energy);
      tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
	  //printf("energy gain=%f\n",tcGain.energy);
      W_tc = exp(-1.0*ffRef.beta*(tcGain.energy + tcLose.energy));
   }
}

inline double MoleculeTransfer::GetCoeff() const
{
#if ENSEMBLE == GEMC
   return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) /
      (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
      boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
   if (sourceBox == mv::BOX0) //Delete case
   {
      return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
         boxDimRef.volInv[sourceBox] *
         exp(-BETA * molRef.kinds[kindIndex].chemPot);
   }
   else //Insertion case
   {
      return boxDimRef.volume[destBox]/
         (double)(molLookRef.NumKindInBox(kindIndex, destBox)+1) *
         exp(BETA * molRef.kinds[kindIndex].chemPot);
   }
#endif
}

inline void MoleculeTransfer::Accept(const uint rejectState, const uint step)
{
   bool result;
   //If we didn't skip the move calculation
   if(rejectState == mv::fail_state::NO_FAIL)
   {
      double molTransCoeff = GetCoeff();
      double Wo = oldMol.GetWeight();
      double Wn = newMol.GetWeight();
      double Wrat = Wn / Wo * W_tc;

      result = prng() < molTransCoeff * Wrat;
      if (result)
      {
         //std::cout << "ACCEPTED\n";
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
         sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
         sysPotRef.boxVirial[sourceBox].tc += tcLose.virial;
         sysPotRef.boxVirial[destBox].tc += tcGain.virial;

         //Add rest of energy.
         sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
         sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
         sysPotRef.boxVirial[sourceBox].inter -= oldVirial;
         sysPotRef.Total();
#ifdef DEBUG_MOVES
         double oldBond = 0.0;
         double oldNonbond = 0.0;
         double newBond = 0.0;
         double newNonbond = 0.0;
         calcEnRef.MoleculeIntra(oldBond, oldNonbond, molIndex, sourceBox);
         double oldInter = calcEnRef.CheckMoleculeInter(molIndex, sourceBox);
#endif

         newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
         comCurrRef.SetNew(molIndex, destBox);
         molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);

         double newVirial = calcEnRef.MoleculeVirial(molIndex, destBox);
         sysPotRef.boxVirial[destBox].inter += newVirial;
         sysPotRef.Total();
         
#ifdef DEBUG_MOVES
	 std::cout << "=============================" 
		   << std::endl << "STEP: " << step << std::endl;
	 SystemPotential tmp = calcEnRef.SystemTotal();
	 for (uint b = 0; b < BOX_TOTAL; ++b)
	 {
	    std::cout << "BOX " << b << ":" << std::endl;
	    for (uint k = 0; k < molRef.kindsCount; ++k)
	    {
	       if (k!=0)
		  std::cout << " ";
	       std::cout << "k" << k << ":" << molLookRef.NumKindInBox(k, b);
	    }
	    std::cout << std::endl;
	    if (b == sourceBox)
	       std::cout << "Delta Individual TC : " << tcLose.energy;
	    else
	       std::cout << "Delta Individual TC : " << tcGain.energy;
	    std::cout << " Current TC : " << tmp.boxEnergy[b].tc
		      << "    " 
		      << "Recalc TC : " << sysPotRef.boxEnergy[b].tc
		      << std::endl;
	 }
	 std::cout  << "=============================" 
		    << std::endl << std::endl;
#endif
      }
   }
   else  //else we didn't even try because we knew it would fail
      result = false;
    subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, sourceBox);
   moveSetRef.Update(result, subPick, step);
}

#endif

#endif

