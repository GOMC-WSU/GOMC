/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef STATIC_VALS_H
#define STATIC_VALS_H

//General includes
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For VARIABLE_<QUANTITY> conditional defines

//Initialization variables


//Member variables
#include "Forcefield.h"
#include "SimEventFrequency.h"
#include "MoleculeLookup.h"
#include "Molecules.h"

class Setup;
class System;

class StaticVals
{
public:
  StaticVals(Setup & set);
  ~StaticVals() {};
  void Init(Setup & set, System& sys);
  void InitOver(Setup & set, System& sys);
  void IsBoxOrthogonal(config_setup::Volume const& vol);
  void IsBoxOrthogonal(const double cellAngle[][3]);

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  double pressure;
  uint kindOfGEMC;
  bool fixVolBox0;
#endif
  bool isOrthogonal;
  bool multiParticleEnabled;

  Forcefield forcefield;
  SimEventFrequency simEventFreq;
  //All the static molecule info --  kind, start index
  Molecules mol;

  double movePerc[mv::MOVE_KINDS_TOTAL];
  double totalPerc;
  config_setup::MEMCVal  intraMemcVal;
  config_setup::FreeEnergy  freeEnVal;

  //Only include these variables if they're static for this ensemble...
#ifndef  VARIABLE_PARTICLE_NUMBER
  MoleculeLookup molLookup;
#endif
#ifdef  VARIABLE_PARTICLE_NUMBER
  config_setup::MEMCVal   memcVal;
  config_setup::CFCMCVal  cfcmcVal;
  config_setup::TargetSwapCollection targetedSwapVal;
#endif

  bool IsEquil(const uint step)
  {
    return step >= simEventFreq.tillEquil;
  }
  bool DoAdjust(const uint move)
  {
    return move % simEventFreq.perAdjust == 0;
  }
  uint GetPerAdjust() const
  {
    return simEventFreq.perAdjust;
  }
  double AcceptPercent(const uint tempAccept)
  {
    return (double)(tempAccept) / (double)(simEventFreq.perAdjust);
  }

  void InitMovePercents(config_setup::MovePercents const& percent);
};

#endif /*STATIC_VALS_H*/
