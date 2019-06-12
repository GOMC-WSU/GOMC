/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SYSTEM_H
#define SYSTEM_H

#include "EnsemblePreprocessor.h" //For VARIABLE_<QUANTITY> conditional defines
#include "CalculateEnergy.h"
#include "EwaldCached.h"
#include "Ewald.h"
#include "NoEwald.h"

//Member variables
#include "EnergyTypes.h"
#include "Coordinates.h"
#include "PRNG.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "MoleculeLookup.h"
#include "MoveSettings.h"
#include "CellList.h"
#include "Clock.h"
#include "CheckpointSetup.h"

//Initialization variables
class Setup;
class StaticVals;
class MoveBase;

class System
{
public:
  explicit System(StaticVals& statics);

  void Init(Setup const& setupData, ulong & startStep);

  //Runs move, picked at random
  void ChooseAndRunMove(const uint step);

  // Recalculate Trajectory
  void RecalculateTrajectory(Setup & set, uint frameNum);

  //print move time
  void PrintTime();

  //print move time
  void PrintAcceptance();

  // return ewald
  Ewald * GetEwald()
  {
    return calcEwald;
  }

#ifdef VARIABLE_VOLUME
  BoxDimensions * BoxDim(const bool isOrthogonal)
  {
    boxDimensions = NULL;
    if(isOrthogonal) {
      boxDimensions = new BoxDimensions();
    } else {
      boxDimensions = new BoxDimensionsNonOrth();
    }
    return boxDimensions;
  }
#endif



  //NOTE:
  //This must also come first... as subsequent values depend on obj.
  //That may be in here, i.e. Box Dimensions
  StaticVals & statV;

  //NOTE:
  //Important! These must come first, as other objects may depend
  //on their val for init!
  //Only include these variables if they vary for this ensemble...
#ifdef VARIABLE_VOLUME
  BoxDimensions *boxDimensions;
#endif
#ifdef  VARIABLE_PARTICLE_NUMBER
  MoleculeLookup molLookup;
#endif

  //Use as we don't know where they are...
  BoxDimensions & boxDimRef;
  MoleculeLookup & molLookupRef;

  MoveSettings moveSettings;
  SystemPotential potential;
  Coordinates coordinates;
  XYZArray atomForceRef;
  XYZArray molForceRef;
  XYZArray atomForceRecRef;
  XYZArray molForceRecRef;
  COM com;

  CalculateEnergy calcEnergy;
  Ewald *calcEwald;
  CellList cellList;
  PRNG prng;

  CheckpointSetup checkpointSet;

  //Procedure to run once move is picked... can also be called directly for
  //debugging...
  void RunMove(uint majKind, real draw, const uint step);

  ~System();


private:
  void InitMoves(Setup const& set);
  void PickMove(uint & kind, real & draw);
  uint SetParams(const uint kind, const real draw);
  uint Transform(const uint kind);
  void CalcEn(const uint kind);
  void Accept(const uint kind, const uint rejectState, const uint step);

  real moveTime[mv::MOVE_KINDS_TOTAL];
  MoveBase * moves[mv::MOVE_KINDS_TOTAL];
  Clock time;
};

#endif /*SYSTEM_H*/
