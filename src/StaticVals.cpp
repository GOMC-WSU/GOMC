/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "StaticVals.h"
#include "ConfigSetup.h" //For types directly read from config. file
#include "GeomLib.h"
#include "Setup.h" //For source of setup data.

using namespace geom;

void StaticVals::Init(Setup & set, System& sys)
{
  //Standard inits
  simEventFreq.Init(set.config.sys.step);
  forcefield.Init(set);
  mol.Init(set, forcefield, sys);
#ifndef VARIABLE_VOLUME
  boxDimensions->Init(set.config.in.restart, set.config.sys.volume,
                      set.pdb.cryst, forcefield);
#endif
#ifndef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(mol, set.pdb.atoms);
#endif
  InitMovePercents(set.config.sys.moves);
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  kindOfGEMC = set.config.sys.gemc.kind;
  pressure = set.config.sys.gemc.pressure;
  fixVolBox0 = set.config.sys.volume.cstVolBox0;
#endif
}

void StaticVals::InitOver(Setup & set, System& sys)
{
  mol.~Molecules();
  mol.Init(set, forcefield, sys);
}

void StaticVals::InitMovePercents(config_setup::MovePercents const& perc)
{
  totalPerc = 0.0;
  for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
    switch(m) {
    case mv::DISPLACE:
      movePerc[m] = perc.displace;
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
#ifdef VARIABLE_VOLUME
    case mv::VOL_TRANSFER :
      movePerc[m] = perc.volume;
      break;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
    case mv::MOL_TRANSFER :
      movePerc[m] = perc.transfer;
      break;
    case mv::MEMC :
      movePerc[m] = perc.memc;
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

void StaticVals::IsBoxOrthogonal(config_setup::Volume const& vol)
{
  real cosAngle[BOX_TOTAL][3];
  real orthogonal[BOX_TOTAL];

  for (uint b = 0; b < BOX_TOTAL; b++) {
    real cellLengthX = vol.axis[b].Length(0);
    real cellLengthY = vol.axis[b].Length(1);
    real cellLengthZ = vol.axis[b].Length(2);
    //Find Cosine Angle of alpha, beta and gamma
    cosAngle[b][0] = Dot(vol.axis[b].Get(1), vol.axis[b].Get(2)) /
                     (cellLengthY * cellLengthZ);
    cosAngle[b][1] = Dot(vol.axis[b].Get(0), vol.axis[b].Get(2)) /
                     (cellLengthX * cellLengthZ);
    cosAngle[b][2] = Dot(vol.axis[b].Get(0), vol.axis[b].Get(1)) /
                     (cellLengthX * cellLengthY);

    orthogonal[b] = ((cosAngle[b][0] == 0.0) &&
                     (cosAngle[b][1] == 0.0) &&
                     (cosAngle[b][2] == 0.0));
    isOrthogonal = (isOrthogonal && orthogonal[b]);
  }
}


StaticVals::StaticVals(Setup & set) : memcVal(set.config.sys.memcVal),
  intraMemcVal(set.config.sys.intraMemcVal)
{
  isOrthogonal = true;
  IsBoxOrthogonal(set.config.sys.volume);
#ifndef VARIABLE_VOLUME
  boxDimensions = NULL;
  if(isOrthogonal) {
    boxDimensions = new BoxDimensions();
  } else {
    boxDimensions = new BoxDimensionsNonOrth();
  }
#endif
}

StaticVals::~StaticVals()
{
#ifndef VARIABLE_VOLUME
  delete boxDimensions;
#endif
}

