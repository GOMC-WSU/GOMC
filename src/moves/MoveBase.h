/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef TRANSFORMABLE_BASE_H
#define TRANSFORMABLE_BASE_H

#include "BasicTypes.h" //For uint.
#include "Molecules.h" //For start
#include "BoxDimensions.h" //For pbc wrapping
#include "BoxDimensionsNonOrth.h"
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

class MoveBase
{
public:

  MoveBase(System & sys, StaticVals const& statV) :
    boxDimRef(sys.boxDimRef), moveSetRef(sys.moveSettings),
    sysPotRef(sys.potential),
    calcEnRef(sys.calcEnergy), comCurrRef(sys.com),
    coordCurrRef(sys.coordinates), prng(sys.prng), molRef(statV.mol),
    BETA(statV.forcefield.beta), ewald(statV.forcefield.ewald),
    cellList(sys.cellList), molRemoved(false),
    atomForceRef(sys.atomForceRef),
    molForceRef(sys.molForceRef),
    atomForceRecRef(sys.atomForceRecRef),
    molForceRecRef(sys.molForceRecRef)
  {
    atomForceNew.Init(sys.atomForceRef.Count());
    molForceNew.Init(sys.molForceRef.Count());
    calcEwald = sys.GetEwald();
    molRemoved = false;
    overlap = false;
    multiParticleEnabled = sys.statV.multiParticleEnabled;
  }

  //Based on the random draw, determine the move kind, box, and
  //(if necessary) molecule kind.
  virtual uint Prep(const real subDraw, const real movPerc) = 0;

  //Note, in general this function is responsible for generating the new
  //configuration to test.
  virtual uint Transform() = 0;

  //In general, this function is responsible for calculating the
  //energy of the system for the new trial configuration.
  virtual void CalcEn() = 0;

  //This function carries out actions based on the internal acceptance state.
  virtual void Accept(const uint rejectState, const uint step) = 0;

  //This function carries out actions based on the internal acceptance state and
  //molecule kind

  //This function print the internal acceptance state for each molecule kind
  virtual void PrintAcceptKind() = 0;

  virtual ~MoveBase() {}

protected:
  uint subPick;
  //If a single molecule move, this is set by the target.
  MoveSettings & moveSetRef;
  SystemPotential & sysPotRef;
  Coordinates & coordCurrRef;
  COM & comCurrRef;
  CalculateEnergy & calcEnRef;
  Ewald * calcEwald;
  XYZArray& atomForceRef;
  XYZArray atomForceNew;
  XYZArray& molForceRef;
  XYZArray molForceNew;
  XYZArray& atomForceRecRef;
  XYZArray& molForceRecRef;

  PRNG & prng;
  BoxDimensions & boxDimRef;
  Molecules const& molRef;
  const real BETA;
  const bool ewald;
  CellList& cellList;
  bool molRemoved, fixBox0, overlap;
  bool multiParticleEnabled;
};

//Data needed for transforming a molecule's position via inter or intrabox
//moves.
class MolTransformBase
{
protected:
  uint GetBoxAndMol(PRNG & prng, Molecules const& molRef,
                    const real subDraw, const real movPerc);
  void ReplaceWith(MolTransformBase const& other);

  //Box, molecule, and molecule kind
  uint b, m, mk;
  uint pStart, pLen;
  //Position
  XYZArray newMolPos;
};

inline uint MolTransformBase::GetBoxAndMol(PRNG & prng, Molecules const& molRef,
    const real subDraw, const real movPerc)
{
#if ENSEMBLE == GCMC
  b = mv::BOX0;
  uint state = prng.PickMol(m, mk, b, subDraw, movPerc);
#else
  uint state = prng.PickMolAndBox(m, mk, b, subDraw, movPerc);
#endif
  pStart = pLen = 0;
  if(state == mv::fail_state::NO_FAIL) {
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

#endif /*TRANSFORMABLE_BASE_H*/
