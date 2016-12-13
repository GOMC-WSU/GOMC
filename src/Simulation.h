/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
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

class Simulation
{
public:
  explicit Simulation(char const*const configFileName);
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
};

#endif /*SIMULATION_H*/
