/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef SIMULATION_H
#define SIMULATION_H
// Member vars
#include "BasicTypes.h"
#include "CPUSide.h"
#include "GOMC_Config.h" //For PT
#include "ParallelTemperingPreprocessor.h"
#include "ParallelTemperingUtilities.h"
#include "StaticVals.h"
#include "System.h"

class Simulation {
public:
  explicit Simulation(char const *const configFileName,
                      MultiSim const *const &multisim = NULL);
  ~Simulation();

  void RunSimulation(void);
  bool RecalculateAndCheck(void);
#if GOMC_GTEST
  ulong GetTrueStep();
  ulong GetRunSteps();
  MoleculeLookup &GetMolLookup();
  MoveSettings &GetMoveSettings();
  Coordinates &GetCoordinates();
  Velocity &GetVelocities();
  ExtendedSystem &GetXSC();
  BoxDimensions &GetBoxDim();
  SystemPotential &GetSystemEnergy(void);
  PRNG &GetPRNG();
  Molecules &GetMolecules();
#endif
private:
  StaticVals *staticValues;
  System *system;
  CPUSide *cpu;
  ulong startStep;
  ulong totalSteps;
  Setup set;
  std::vector<ulong> frameSteps;
  uint remarksCount;
  MultiSim const *const &ms;
#if GOMC_LIB_MPI
  ParallelTemperingUtilities *PTUtils;
  std::vector<bool> exchangeResults;
  int parity;
#endif
};

#endif /*SIMULATION_H*/
