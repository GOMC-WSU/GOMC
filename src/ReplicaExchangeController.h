/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef REPLICAEXCHANGECONTROLLER_H
#define REPLICAEXCHANGECONTROLLER_H

#include <stdlib.h>
#include "Simulation.h"
#include "MersenneTwister.h"

using namespace std; 

class ReplicaExchangeController
{
public:
    explicit ReplicaExchangeController(vector<Simulation*>*);
    ~ReplicaExchangeController();
    void runMultiSim();
    double calcDelta(int j);
    void exchange(int j);

private:
    vector<Simulation*>* simsRef;
    ulong exchangeRate;
    ulong totalSteps;
    ulong roundedUpDivison;
    double swapperForT_in_K;
    double swapperForBeta;
    CPUSide * swapperForCPUSide;
    int parityOfSwaps;
    double checkerForIncreasingMontonicityOfTemp;
    MTRand rand;
};

#endif