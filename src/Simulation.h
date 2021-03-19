/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SIMULATION_H
#define SIMULATION_H

//Member vars
#include "CPUSide.h"
#include "System.h"
#include "StaticVals.h"
#include "BasicTypes.h"
#include "GOMC_Config.h"    //For PT
#include "ParallelTemperingPreprocessor.h"
#include "ParallelTemperingUtilities.h"

class Simulation
{
public:
  explicit Simulation(char const*const configFileName, MultiSim const*const& multisim = NULL);
  ~Simulation();

  void RunSimulation(void);
  bool RecalculateAndCheck(void);

  #if GOMC_GTEST_MPI
      double GetSystemEnergy(void);
      void ExchangeReplicas(int worldRank);
      Coordinates& getCoordinates(void);
      COM& getCOMs(void);

  #endif

private:
  StaticVals * staticValues;
  System * system;
  CPUSide * cpu;
  ulong totalSteps;
  Setup set;
  std::vector<ulong> frameSteps;
  uint remarksCount;
  ulong startStep;
  MultiSim const*const& ms;
#if GOMC_LIB_MPI
  ParallelTemperingUtilities * PTUtils;
  std::vector<bool> exchangeResults;
  int parity;
#endif
};

#endif /*SIMULATION_H*/
