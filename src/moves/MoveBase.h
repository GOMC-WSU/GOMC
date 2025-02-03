/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef TRANSFORMABLE_BASE_H
#define TRANSFORMABLE_BASE_H

#include "BasicTypes.h"    //For uint.
#include "BoxDimensions.h" //For pbc wrapping
#include "BoxDimensionsNonOrth.h"
#include "COM.h"
#include "CalculateEnergy.h"
#include "Coordinates.h"
#include "EnergyTypes.h"
#include "Ewald.h"
#include "EwaldCached.h"
#include "Forcefield.h"
#include "GOMCEventsProfile.h" // for NVTX profiling
#include "MolPick.h"
#include "Molecules.h" //For start
#include "MoveConst.h"
#include "MoveSettings.h"
#include "NoEwald.h"
#include "StaticVals.h"
#include "System.h"
#include "Velocity.h"
#include "XYZArray.h" //Parent class

class MoveBase {
public:
  MoveBase(System &sys, StaticVals const &statV)
      : moveSetRef(sys.moveSettings), sysPotRef(sys.potential),
        coordCurrRef(sys.coordinates), comCurrRef(sys.com),
        calcEnRef(sys.calcEnergy), atomForceRef(sys.atomForceRef),
        molForceRef(sys.molForceRef), atomForceRecRef(sys.atomForceRecRef),
        molForceRecRef(sys.molForceRecRef), velocity(sys.vel), prng(sys.prng),
        boxDimRef(sys.boxDimRef), molRef(statV.mol),
        BETA(statV.forcefield.beta), ewald(statV.forcefield.ewald),
        cellList(sys.cellList) {
    calcEwald = sys.GetEwald();
    molRemoved = false;
    overlap = false;
    multiParticleEnabled = sys.statV.multiParticleEnabled;
  }

  // Based on the random draw, determine the move kind, box, and
  //(if necessary) molecule kind.
  virtual uint Prep(const double subDraw, const double movPerc) = 0;

  // Setup the picked box and molecule to perform relaxation in NeMTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0) = 0;

  // Note, in general this function is responsible for generating the new
  // configuration to test.
  virtual uint Transform() = 0;

  // In general, this function is responsible for calculating the
  // energy of the system for the new trial configuration.
  virtual void CalcEn() = 0;

  // This function carries out actions based on the internal acceptance state.
  virtual void Accept(const uint rejectState, const ulong step) = 0;

  // This function carries out actions based on the internal acceptance state
  // and molecule kind

  // This function print the internal acceptance state for each molecule kind
  virtual void PrintAcceptKind() = 0;

  virtual ~MoveBase() {}

protected:
  uint subPick;
  // If a single molecule move, this is set by the target.
  MoveSettings &moveSetRef;
  SystemPotential &sysPotRef;
  Coordinates &coordCurrRef;
  COM &comCurrRef;
  CalculateEnergy &calcEnRef;
  Ewald *calcEwald;
  XYZArray &atomForceRef;
  XYZArray &molForceRef;
  XYZArray &atomForceRecRef;
  XYZArray &molForceRecRef;
  Velocity &velocity;

  PRNG &prng;
  BoxDimensions &boxDimRef;
  Molecules const &molRef;
  const double BETA;
  const bool ewald;
  CellList &cellList;
  bool molRemoved, fixBox0, overlap;
  bool multiParticleEnabled;
};

// Data needed for transforming a molecule's position via inter or intra box
// moves.
class MolTransformBase {
protected:
  uint GetBoxAndMol(PRNG &prng, Molecules const &molRef, const double subDraw,
                    const double movPerc);
  void ReplaceWith(MolTransformBase const &other);

  // Box, molecule, and molecule kind
  uint b, m, mk;
  uint pStart, pLen;
  // Position
  XYZArray newMolPos;
};

inline uint MolTransformBase::GetBoxAndMol(PRNG &prng, Molecules const &molRef,
                                           const double subDraw,
                                           const double movPerc) {
#if ENSEMBLE == GCMC
  b = mv::BOX0;
  uint state = prng.PickMol(m, mk, b, subDraw, movPerc);
#else
  uint state = prng.PickMolAndBox(m, mk, b, subDraw, movPerc);
#endif
  pStart = pLen = 0;
  if (state == mv::fail_state::NO_FAIL) {
    molRef.GetRangeStartLength(pStart, pLen, m);
    newMolPos.Uninit();
    newMolPos.Init(pLen);
  }
  return state;
}

inline void MolTransformBase::ReplaceWith(MolTransformBase const &other) {
  m = other.m;
  mk = other.mk;
  b = other.b;
  pStart = other.pStart;
  pLen = other.pLen;
  newMolPos = other.newMolPos;
}

#endif /*TRANSFORMABLE_BASE_H*/
