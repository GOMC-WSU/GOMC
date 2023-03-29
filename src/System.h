/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef SYSTEM_H
#define SYSTEM_H

#include "CalculateEnergy.h"
#include "EnsemblePreprocessor.h" //For VARIABLE_<QUANTITY> conditional defines
#include "Ewald.h"
#include "EwaldCached.h"
#include "NoEwald.h"

// Member variables
#include "../lib/Lambda.h"
#include "../lib/StrLib.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "CellList.h"
#include "CheckpointSetup.h"
#include "Clock.h"
#include "Coordinates.h"
#include "EnergyTypes.h"
#include "ExtendedSystem.h"
#include "MoleculeLookup.h"
#include "MoveSettings.h"
#include "PRNG.h"
#include "Random123Wrapper.h"
#include "Velocity.h"

// Initialization variables
class Setup;
class StaticVals;
class MoveBase;
class Lambda;

class System {
public:
  explicit System(StaticVals &statics, Setup &set, ulong &startStep,
                  MultiSim const *const &multisim = NULL);

  void Init(Setup &setupData);

  /* To reinit the checkpointed original molecule starts */
  void InitOver(Setup &set, Molecules &molRef);

  // Runs move, picked at random
  void ChooseAndRunMove(const ulong step);

  // Recalculate Trajectory
  void RecalculateTrajectory(Setup &set, uint frameNum);

  // print move time
  void PrintTime();

  // print move time
  void PrintAcceptance();

  // return ewald
  Ewald *GetEwald() { return calcEwald; }

  // Return the pointer to specific move
  MoveBase *GetMoveObject(const uint moveIndex) { return moves[moveIndex]; }

  BoxDimensions *BoxDim(const bool isOrthogonal) {
    boxDimensions = NULL;
    if (isOrthogonal) {
      boxDimensions = new BoxDimensions();
    } else {
      boxDimensions = new BoxDimensionsNonOrth();
    }
    return boxDimensions;
  }

  // NOTE:
  // This must also come first... as subsequent values depend on obj.
  // That may be in here, i.e. Box Dimensions
  StaticVals &statV;

  // NOTE:
  // Important! These must come first, as other objects may depend
  // on their val for init!
  // Only include these variables if they vary for this ensemble...
  BoxDimensions *boxDimensions;
#ifdef VARIABLE_PARTICLE_NUMBER
  MoleculeLookup molLookup;
#endif

  // Use as we don't know where they are...
  BoxDimensions &boxDimRef;
  MoleculeLookup &molLookupRef;

  PRNG prng;
  Random123Wrapper r123wrapper;

#if GOMC_LIB_MPI
  MultiSim const *const &ms;
  PRNG *prngParallelTemp;
#endif

  MoveSettings moveSettings;
  CellList cellList;
  SystemPotential potential;
  Coordinates coordinates;
  XYZArray atomForceRef;
  XYZArray molForceRef;
  XYZArray atomForceRecRef;
  XYZArray molForceRecRef;
  Lambda lambdaRef;
  COM com;
  ExtendedSystem xsc;
  Velocity vel;

  CalculateEnergy calcEnergy;
  Ewald *calcEwald;

  /* For checkpoint restoration */
  CheckpointSetup checkpointSet;
  bool restartFromCheckpoint;
  ulong &startStepRef, trueStep;
  // Procedure to run once move is picked... can also be called directly for
  // debugging...
  void RunMove(uint majKind, double draw, const ulong step);

  ~System();

private:
  void InitLambda();
  void InitMoves(Setup const &set);
  void PickMove(uint &kind, double &draw);
  uint SetParams(const uint kind, const double draw);
  uint Transform(const uint kind);
  void CalcEn(const uint kind);
  void Accept(const uint kind, const uint rejectState, const ulong step);

  double moveTime[mv::MOVE_KINDS_TOTAL];
  MoveBase *moves[mv::MOVE_KINDS_TOTAL];
  Clock time;
};

#endif /*SYSTEM_H*/
