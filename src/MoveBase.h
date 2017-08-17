/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.0
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef TRANSFORMABLE_BASE_H
#define TRANSFORMABLE_BASE_H

#include "BasicTypes.h" //For uint.
#include "Molecules.h" //For start
#include "BoxDimensions.h" //For pbc wrapping
#include "XYZArray.h" //Parent class
#include "MoveSettings.h"
#include "Coordinates.h"
#include "EnergyTypes.h"
#include "COM.h"
#include "MoveConst.h"
#include "System.h"
#include "StaticVals.h"
#include "CalculateEnergy.h"
#include "EwaldCached.h"
#include "Ewald.h"
#include "NoEwald.h"
#include "MolPick.h"
#include "Forcefield.h"
#ifdef _OPENMP
#include <omp.h>
#endif

class MoveBase
{
 public:

   MoveBase(System & sys, StaticVals const& statV) :
   boxDimRef(sys.boxDimRef), moveSetRef(sys.moveSettings),
     sysPotRef(sys.potential),
     calcEnRef(sys.calcEnergy), comCurrRef(sys.com),
     coordCurrRef(sys.coordinates), prng(sys.prng), molRef(statV.mol),
     BETA(statV.forcefield.beta), ewald(statV.forcefield.ewald),
     cellList(sys.cellList), molRemoved(false)
   {
      calcEwald = sys.GetEwald();
#if ENSEMBLE == GEMC || ENSEMBLE == NPT     
      fixBox0 = statV.fixVolBox0;
#endif
   }

    //Based on the random draw, determine the move kind, box, and
    //(if necessary) molecule kind.
    virtual uint Prep(const double subDraw, const double movPerc) = 0;

    //Note, in general this function is responsible for generating the new
    //configuration to test.
    virtual uint Transform() = 0;

    //In general, this function is responsible for calculating the
    //energy of the system for the new trial configuration.
    virtual void CalcEn() = 0;

    //This function carries out actions based on the internal acceptance state.
    virtual void Accept(const uint rejectState, const uint step) = 0;

 protected:
    uint subPick;
    //If a single molecule move, this is set by the target.
    MoveSettings & moveSetRef;
    SystemPotential & sysPotRef;
    Coordinates & coordCurrRef;
    COM & comCurrRef;
    CalculateEnergy & calcEnRef;
    EwaldCached * calcEwald;

    PRNG & prng;
    BoxDimensions & boxDimRef;
    Molecules const& molRef;
    const double BETA;
    const bool ewald;
    CellList& cellList;
    bool molRemoved, fixBox0;
};

//Data needed for transforming a molecule's position via inter or intrabox
//moves.
class MolTransformBase
{
 protected:
   uint GetBoxAndMol(PRNG & prng, Molecules const& molRef,
		     const double subDraw, const double movPerc);
   void ReplaceWith(MolTransformBase const& other);

   //Box, molecule, and molecule kind
   uint b, m, mk;
   uint pStart, pLen;
   //Position
   XYZArray newMolPos;
};

inline uint MolTransformBase::GetBoxAndMol
(PRNG & prng, Molecules const& molRef,
 const double subDraw, const double movPerc)
{
#if ENSEMBLE == GCMC
   b = mv::BOX0;
   uint state = prng.PickMol(m, mk, b, subDraw, movPerc);
#else
   uint state = prng.PickMolAndBox(m, mk, b, subDraw, movPerc);
#endif
   pStart = pLen = 0;
   if(state == mv::fail_state::NO_FAIL)
   {
     molRef.GetRangeStartLength(pStart, pLen, m);
     newMolPos.Uninit();
     newMolPos.Init(pLen);
   }
   return state;
}

inline void MolTransformBase::ReplaceWith(MolTransformBase const& other)
{
   m = other.m;
   mk = other.mk;
   b = other.b;
   pStart = other.pStart;
   pLen = other.pLen;
   newMolPos = other.newMolPos;
}

class Rotate;

class Translate : public MoveBase, public MolTransformBase
{
 public:

   Translate(System &sys, StaticVals const& statV) : MoveBase(sys, statV) {}

   virtual uint Prep(const double subDraw, const double movPerc);
   uint ReplaceRot(Rotate const& other);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint rejectState, const uint step);
 private:
   Intermolecular inter_LJ, inter_Real, recip;
   XYZ newCOM;
};

inline uint Translate::Prep(const double subDraw, const double movPerc)
{ return GetBoxAndMol(prng, molRef, subDraw, movPerc); }

inline uint Translate::Transform()
{
   subPick = mv::GetMoveSubIndex(mv::DISPLACE, b);
   coordCurrRef.TranslateRand(newMolPos, newCOM, pStart, pLen,
			      m, b, moveSetRef.Scale(subPick));
   return mv::fail_state::NO_FAIL;
}

inline void Translate::CalcEn()
{
   cellList.RemoveMol(m, b, coordCurrRef);
   molRemoved = true;

   //calculate LJ interaction and real term of electrostatic interaction
   calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
   //calculate reciprocate term of electrostatic interaction
   recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);

}

inline void Translate::Accept(const uint rejectState, const uint step)
{
   bool res =false;
   if (rejectState == mv::fail_state::NO_FAIL)
   {
      double pr = prng();
      res = pr < exp(-BETA * (inter_LJ.energy + inter_Real.energy +
			      recip.energy));
   }
   bool result = (rejectState == mv::fail_state::NO_FAIL) && res;

   if (result)
   {
      //Set new energy.
      // setting energy and virial of LJ interaction
      sysPotRef.boxEnergy[b].inter += inter_LJ.energy;
      // setting energy and virial of coulomb interaction
      sysPotRef.boxEnergy[b].real += inter_Real.energy;
      // setting energy and virial of recip term
      sysPotRef.boxEnergy[b].recip += recip.energy;;

      //Copy coords
      newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
      comCurrRef.Set(m, newCOM);
      calcEwald->UpdateRecip(b);

      sysPotRef.Total();
   }

   if(molRemoved)
   {
     // It means that Recip energy is calculated and move not accepted
     if(!result)
     {
       calcEwald->RestoreMol(m);
     }
     cellList.AddMol(m, b, coordCurrRef);
     molRemoved = false;
   }

   subPick = mv::GetMoveSubIndex(mv::DISPLACE, b);
   moveSetRef.Update(result, subPick, step);
}

class Rotate : public MoveBase, public MolTransformBase
{
 public:
   Rotate(System &sys, StaticVals const& statV) : MoveBase(sys, statV) {}

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual uint Transform();
   virtual void CalcEn();
   virtual void Accept(const uint earlyReject, const uint step);
 private:
   Intermolecular inter_LJ, inter_Real, recip;
};

inline uint Rotate::Prep(const double subDraw, const double movPerc)
{
   uint state = GetBoxAndMol(prng, molRef, subDraw, movPerc);
   if (state == mv::fail_state::NO_FAIL && molRef.NumAtoms(mk)  <= 1)
	 state = mv::fail_state::ROTATE_ON_SINGLE_ATOM;
   return state;
}

inline uint Translate::ReplaceRot(Rotate const& other)
{
   ReplaceWith(other);
   return mv::fail_state::NO_FAIL;
}

inline uint Rotate::Transform()
{
   subPick = mv::GetMoveSubIndex(mv::ROTATE, b);
   coordCurrRef.RotateRand(newMolPos, pStart, pLen, m, b,
			   moveSetRef.Scale(subPick));
   return mv::fail_state::NO_FAIL;
}

inline void Rotate::CalcEn()
{
   cellList.RemoveMol(m, b, coordCurrRef);
   molRemoved = true;

   //calculate LJ interaction and real term of electrostatic interaction
   calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
   //calculate reciprocate term of electrostatic interaction
   recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
}

inline void Rotate::Accept(const uint rejectState, const uint step)
{
   bool res =false;

   if (rejectState == mv::fail_state::NO_FAIL)
   {
      double pr = prng();
      res = pr < exp(-BETA * (inter_LJ.energy + inter_Real.energy +
			      recip.energy));
   }
   bool result = (rejectState == mv::fail_state::NO_FAIL) && res;

   if (result)
   {
      //Set new energy.
      // setting energy and virial of LJ interaction
      sysPotRef.boxEnergy[b].inter += inter_LJ.energy;
      // setting energy and virial of coulomb interaction
      sysPotRef.boxEnergy[b].real += inter_Real.energy;
      // setting energy and virial of recip term
      sysPotRef.boxEnergy[b].recip += recip.energy;

      //Copy coords
      newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
      calcEwald->UpdateRecip(b);

      sysPotRef.Total();
   }

   if(molRemoved)
   {
     // It means that Recip energy is calculated and move not accepted
     if(!result)
     {
       calcEwald->RestoreMol(m);
     }

     cellList.AddMol(m, b, coordCurrRef);
     molRemoved = false;
   }

   subPick = mv::GetMoveSubIndex(mv::ROTATE, b);
   moveSetRef.Update(result, subPick, step);
}

#if ENSEMBLE == GEMC || ENSEMBLE == NPT

class VolumeTransfer : public MoveBase
{
 public:
   VolumeTransfer(System &sys, StaticVals const& statV);

   virtual uint Prep(const double subDraw, const double movPerc);
   virtual void CalcEn();
   virtual uint Transform();
   double GetCoeff() const;
   virtual void Accept(const uint rejectState, const uint step);
 private:
   //Note: This is only used for GEMC-NPT
   uint bPick[BOX_TOTAL], subPick, subPickT[BOX_TOTAL];
   SystemPotential sysPotNew;
   BoxDimensions newDim;
   Coordinates newMolsPos;
   COM newCOMs;
   MoleculeLookup & molLookRef;
   const uint GEMC_KIND;
   const double PRESSURE;
   bool regrewGrid;
};

inline VolumeTransfer::VolumeTransfer(System &sys, StaticVals const& statV)  :
		      MoveBase(sys, statV), molLookRef(sys.molLookupRef),
		      newDim(sys.boxDimRef), newMolsPos(boxDimRef, newCOMs,
							sys.molLookupRef,
							sys.prng, statV.mol),
		      newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef,
			      statV.mol), GEMC_KIND(statV.kindOfGEMC),
		      PRESSURE(statV.pressure), regrewGrid(false)
{
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(statV.mol.count);
}

inline uint VolumeTransfer::Prep(const double subDraw, const double movPerc)
{
   uint state = mv::fail_state::NO_FAIL;
   if (GEMC_KIND == mv::GEMC_NVT)
   {
      subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER);
      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	bPick[b] = b;
      }
   }
   if (GEMC_KIND == mv::GEMC_NPT)
   {
#if ENSEMBLE == NPT
      prng.PickBox(bPick[0], subDraw, movPerc);
#else
      prng.PickBoxPair(bPick[0], bPick[1], subDraw, movPerc);
#endif
      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	subPickT[bPick[b]] = mv::GetMoveSubIndex(mv::VOL_TRANSFER, bPick[b]);
      }
   }
   newDim = boxDimRef;
   coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
   comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
   return state;
}

inline uint VolumeTransfer::Transform()
{
   uint state = mv::fail_state::NO_FAIL;
   //Reinit, if necessary.
   if (GEMC_KIND == mv::GEMC_NVT)
   {
      double max = moveSetRef.Scale(subPick);
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs, newDim,
					   comCurrRef, max);
   }
   else
   {
      XYZ scale[BOX_TOTAL];
      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	if (state == mv::fail_state::NO_FAIL)
	{
	  if ((bPick[b] == 0) && fixBox0)
	     continue;

	   double max = moveSetRef.Scale(subPickT[bPick[b]]);
	   double delta = prng.Sym(max);
	   state =  boxDimRef.ShiftVolume(newDim, scale[bPick[b]],
					  bPick[b], delta);
	}
      }

      if (state == mv::fail_state::NO_FAIL)
      {
	 for (uint b = 0; b < BOX_TOTAL; b++)
	 {
	    if ((bPick[b] == 0) && fixBox0)
	       continue;

	    coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef, 
				      newDim, bPick[b], scale[bPick[b]]);
	 }
      }
   }
   return state;
}

inline void VolumeTransfer::CalcEn()
{
   cellList.GridAll(newDim, newMolsPos, molLookRef);
   regrewGrid = true;

    //back up cached fourier term
   calcEwald->exgMolCache();
   sysPotNew = sysPotRef;
   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      if ((bPick[b] == 0) && fixBox0)
	continue;

      //calculate new K vectors
      calcEwald->RecipInit(bPick[b], newDim);
      //setup reciprocate terms
      calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
      //calculate LJ interaction and real term of electrostatic interaction
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newCOMs, newDim,
				     bPick[b]);
      //calculate reciprocate term of electrostatic interaction
      sysPotNew.boxEnergy[bPick[b]].recip = calcEwald->BoxReciprocal(bPick[b]);
   }
   sysPotNew.Total();  
}



inline double VolumeTransfer::GetCoeff() const
{
   ////Log-volume style shift -- is turned off, at present.
   //
   //return pow(newDim.volume[b_i]/boxDimRef.volume[b_i],
   //	      (double)molLookRef.NumInBox(b_i)+1) *
   //  pow(newDim.volume[b_ii]/boxDimRef.volume[b_ii],
   //	 (double)molLookRef.NumInBox(b_ii)+1);
   double coeff = 1.0;
   if (GEMC_KIND == mv::GEMC_NVT)
   {
      for (uint b = 0; b < BOX_TOTAL; ++b)
      {
	 coeff *= pow(newDim.volume[b]/boxDimRef.volume[b],
		      (double)molLookRef.NumInBox(b));
      }
   }
   else
   {
      for (uint b = 0; b < BOX_TOTAL; ++b)
      {
	if ((bPick[b] == 0) && fixBox0)
	  continue;

	coeff *= pow(newDim.volume[b]/boxDimRef.volume[b],
		     (double)molLookRef.NumInBox(b)) *
	  exp(-BETA * PRESSURE * (newDim.volume[b]-boxDimRef.volume[b]));
      }

   }
   return coeff;
}

inline void VolumeTransfer::Accept(const uint rejectState, const uint step)
{
   double volTransCoeff = GetCoeff();
   double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
   double accept = volTransCoeff * uBoltz;
   bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
   if (result)
   {
      //Set new energy.
      sysPotRef = sysPotNew;
      //Swap... next time we'll use the current members.
      //NOTE:
      //This will be less efficient for NPT, but necessary evil.
      swap(coordCurrRef, newMolsPos);
      swap(comCurrRef, newCOMs);
      boxDimRef = newDim;

      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	 if ((bPick[b] == 0) && fixBox0)
	   continue;
	 
	 calcEwald->UpdateRecip(b);
	 calcEwald->UpdateRecipVec(b);
      }
   }
   else if (rejectState == mv::fail_state::NO_FAIL && regrewGrid)
   {
      cellList.GridAll(boxDimRef, coordCurrRef, molLookRef);
      regrewGrid = false;

      calcEwald->exgMolCache();
   }

   if (GEMC_KIND == mv::GEMC_NVT)
   {
      subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER, 0);
      moveSetRef.Update(result, subPick, step);
      subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER, 1);
      moveSetRef.Update(result, subPick, step);
   }
   if (GEMC_KIND == mv::GEMC_NPT)
   {
     for (uint b = 0; b < BOX_TOTAL; b++)
     {
       if ((bPick[b] == 0) && fixBox0)
	 continue;

       subPickT[bPick[b]] = mv::GetMoveSubIndex(mv::VOL_TRANSFER, bPick[b]);
       moveSetRef.Update(result, subPickT[bPick[b]], step);
     }
   }

}

#endif

#endif /*TRANSFORMABLE_BASE_H*/
