/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
#if GOMC_LIB_MPI
#include "ParallelTemperingPreprocessor.h"
#endif
class Simulation
{
public:
  explicit Simulation(char const*const configFileName);
  #if GOMC_LIB_MPI
  explicit Simulation(char const*const configFileName, MultiSim *& multisim);
  #endif
  ~Simulation();

  void RunSimulation(void);

#ifndef NDEBUG
  void RunningCheck(const uint step);
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
  ulong step;
};

#endif /*SIMULATION_H*/
