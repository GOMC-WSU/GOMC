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

   double W_tc, W_recip, oldVirial_LJ, oldVirial_Real;
   cbmc::TrialMol oldMol, newMol;
   Intermolecular tcLose, tcGain, recipLose, recipGain;
   MoleculeLookup & molLookRef;
   Forcefield const& ffRef;
};

inline uint MoleculeTransfer::GetBoxPairAndMol
(const double subDraw, const double movPerc)
{
   uint state = prng.PickMolAndBoxPair(molIndex, kindIndex,
				       sourceBox, destBox,
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
   newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef,
			   destBox);
   oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef,
			   sourceBox);
   oldMol.SetCoords(coordCurrRef, pStart);
   W_tc = 1.0;
   return state;
}


inline uint MoleculeTransfer::Transform()
{
   oldVirial_LJ = 0.0; 
   oldVirial_Real = 0.0;
   calcEnRef.MoleculeVirial(oldVirial_LJ, oldVirial_Real, molIndex, sourceBox);
   cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
   subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, sourceBox);
   molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);
   return mv::fail_state::NO_FAIL;
}

inline void MoleculeTransfer::CalcEn()
{
   if (ffRef.useLRC)
   {
      tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
      tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
      W_tc = exp(-1.0*ffRef.beta*(tcGain.energy + tcLose.energy));
      W_recip = 1.0;
      if (newMol.GetWeight() != 0.0)
      {
	 recipGain.energy =
	   calcEwald->SwapDestRecip(newMol, destBox, sourceBox, molIndex);
	 recipLose.energy =
	   calcEwald->SwapSourceRecip(oldMol, sourceBox, molIndex);
	 W_recip = exp(-1.0 * ffRef.beta * (recipGain.energy +
					    recipLose.energy));
      }
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
      double Wrat = Wn / Wo * W_tc * W_recip;

      result = prng() < molTransCoeff * Wrat;
      if (result)
      {
         //Add tail corrections
         sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
         sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
         sysPotRef.boxVirial[sourceBox].tc += tcLose.virial;
         sysPotRef.boxVirial[destBox].tc += tcGain.virial;

         //Add rest of energy.
         sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
         sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
         sysPotRef.boxVirial[sourceBox].inter -= oldVirial_LJ;
	 sysPotRef.boxVirial[sourceBox].real -= oldVirial_Real;
	 sysPotRef.boxEnergy[sourceBox].recip += recipLose.energy;
	 sysPotRef.boxEnergy[destBox].recip += recipGain.energy;

	 //Set coordinates, new COM; shift index to new box's list
         newMol.GetCoords().CopyRange(coordCurrRef, 0, pStart, pLen);
         comCurrRef.SetNew(molIndex, destBox);
         molLookRef.ShiftMolBox(molIndex, sourceBox, destBox,
				kindIndex);
	 cellList.AddMol(molIndex, destBox, coordCurrRef);

	 //Calculate the fresh virial.
         double newVirial_LJ = 0.0, newVirial_Real = 0.0;
	 calcEnRef.MoleculeVirial(newVirial_LJ, newVirial_Real,
				  molIndex,destBox);
         sysPotRef.boxVirial[destBox].inter += newVirial_LJ;
	 sysPotRef.boxVirial[destBox].real += newVirial_Real;

	 //Zero out box energies to prevent small number 
	 //errors in double.
	 if (molLookRef.NumInBox(sourceBox) == 0)
	 {
	    sysPotRef.boxEnergy[sourceBox].Zero();
	    sysPotRef.boxVirial[sourceBox].Zero();
	 }
	 else if (molLookRef.NumInBox(sourceBox) == 1)
	 {
	    sysPotRef.boxEnergy[sourceBox].inter = 0;
	    sysPotRef.boxVirial[sourceBox].inter = 0;
	    sysPotRef.boxEnergy[sourceBox].real = 0;
	    sysPotRef.boxVirial[sourceBox].real = 0;
	 }

	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    calcEwald->UpdateRecip(b);
	 }

	 //Retotal
         sysPotRef.Total();
      }
      else
      {
	 cellList.AddMol(molIndex, sourceBox, coordCurrRef);
	 //when weight is 0, MolDestSwap() will not be executed, thus cos/sin
	 //molRef will not be changed. Also since no memcpy, doing restore
	 //results in memory overwrite
	 if (newMol.GetWeight() != 0.0)
	    calcEwald->RestoreMol(molIndex);
	 
      }

   }
   else  //else we didn't even try because we knew it would fail
      result = false;

   subPick = mv::GetMoveSubIndex(mv::MOL_TRANSFER, sourceBox);
   moveSetRef.Update(result, subPick, step);
}

#endif

#endif
