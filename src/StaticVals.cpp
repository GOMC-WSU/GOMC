/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "StaticVals.h"

#include "ConfigSetup.h" //For types directly read from config. file
#include "GeomLib.h"
#include "Setup.h" //For source of setup data.

using namespace geom;

void StaticVals::Init(Setup &set, System &sys) {
  // Standard inits
  simEventFreq.Init(set.config.sys.step);
  forcefield.Init(set);
  mol.Init(set, forcefield, sys);
#ifndef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(mol, set.pdb.atoms, forcefield,
                 set.config.in.restart.restartFromCheckpoint);
#endif
  InitMovePercents(set.config.sys.moves);
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  kindOfGEMC = set.config.sys.gemc.kind;
  pressure = set.config.sys.gemc.pressure;
  fixVolBox0 = set.config.sys.volume.cstVolBox0;
#endif
}

void StaticVals::InitOver(Setup &set, System &sys) {
  mol.~Molecules();
  mol.Init(set, forcefield, sys);
}

void StaticVals::InitMovePercents(config_setup::MovePercents const &perc) {
  totalPerc = 0.0;
  for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
    switch (m) {
    case mv::DISPLACE:
      movePerc[m] = perc.displace;
      break;
    case mv::MULTIPARTICLE:
      movePerc[m] = perc.multiParticle;
      break;
    case mv::MULTIPARTICLE_BM:
      movePerc[m] = perc.multiParticleBrownian;
      break;
    case mv::ROTATE:
      movePerc[m] = perc.rotate;
      break;
    case mv::INTRA_SWAP:
      movePerc[m] = perc.intraSwap;
      break;
    case mv::REGROWTH:
      movePerc[m] = perc.regrowth;
      break;
    case mv::INTRA_MEMC:
      movePerc[m] = perc.intraMemc;
      break;
    case mv::CRANKSHAFT:
      movePerc[m] = perc.crankShaft;
      break;
    case mv::INTRA_TARGETED_SWAP:
      movePerc[m] = perc.intraTargetedSwap;
      break;
#ifdef VARIABLE_VOLUME
    case mv::VOL_TRANSFER:
      movePerc[m] = perc.volume;
      break;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
    case mv::MOL_TRANSFER:
      movePerc[m] = perc.transfer;
      break;
    case mv::MEMC:
      movePerc[m] = perc.memc;
      break;
    case mv::NE_MTMC:
      movePerc[m] = perc.neMolTransfer;
      break;
    case mv::TARGETED_SWAP:
      movePerc[m] = perc.targetedSwap;
      break;
#endif
#endif
    default:
      movePerc[m] = 0.0;
      break;
    }
    totalPerc += movePerc[m];
  }
  for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++)
    movePerc[m] /= totalPerc;
  totalPerc = 1.0;
}

void StaticVals::IsBoxOrthogonal(config_setup::Volume const &vol) {
  double cosAngle[BOX_TOTAL][3];
  bool orthogonal[BOX_TOTAL];

  for (uint b = 0; b < BOX_TOTAL; b++) {
    double cellLengthX = vol.axis[b].Length(0);
    double cellLengthY = vol.axis[b].Length(1);
    double cellLengthZ = vol.axis[b].Length(2);
    // Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = Dot(vol.axis[b].Get(1), vol.axis[b].Get(2)) /
                     (cellLengthY * cellLengthZ);
    cosAngle[b][1] = Dot(vol.axis[b].Get(0), vol.axis[b].Get(2)) /
                     (cellLengthX * cellLengthZ);
    cosAngle[b][2] = Dot(vol.axis[b].Get(0), vol.axis[b].Get(1)) /
                     (cellLengthX * cellLengthY);

    orthogonal[b] = ((cosAngle[b][0] == 0.0) && (cosAngle[b][1] == 0.0) &&
                     (cosAngle[b][2] == 0.0));
    isOrthogonal &= orthogonal[b];
  }
}

void StaticVals::IsBoxOrthogonal(const double cellAngle[][3]) {
  for (uint b = 0; b < BOX_TOTAL; b++) {
    bool orthogonal = ((cellAngle[b][0] == 90.0) && (cellAngle[b][1] == 90.0) &&
                       (cellAngle[b][2] == 90.0));
    isOrthogonal &= orthogonal;
  }
}

StaticVals::StaticVals(Setup &set)
    : intraMemcVal(set.config.sys.intraMemcVal),
      freeEnVal(set.config.sys.freeEn), memcVal(set.config.sys.memcVal),
      neMTMCVal(set.config.sys.neMTMCVal),
      targetedSwapVal(set.config.sys.targetedSwapCollection),
      intraTargetedSwapVal(set.config.sys.intraTargetedSwapCollection)

{
  multiParticleEnabled = set.config.sys.moves.multiParticleEnabled;
  multiParticleLiquid = set.config.sys.moves.multiParticleLiquid;
  multiParticleGas = set.config.sys.moves.multiParticleGas;
  isOrthogonal = true;
  if (set.config.in.restart.enable) {
    IsBoxOrthogonal(set.pdb.cryst.cellAngle);
  } else {
    IsBoxOrthogonal(set.config.sys.volume);
  }
}
