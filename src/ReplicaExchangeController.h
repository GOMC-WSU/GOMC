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
#include "Writer.h"
#include "Clock.h"

using namespace std; 

class ReplicaExchangeController
{
public:
    explicit ReplicaExchangeController(vector<Simulation*>*);
    ~ReplicaExchangeController();
    void runMultiSim();
    double calc_delta(FILE * fplog, int a, int b, int ap, int bp);
    void exchange(int a, int b);
    void swapSimulations(int a, int b);
    void InitRecordKeeper();
    void print_ind(FILE * fplog, const char *leg, int n, int *ind, bool *bEx);
    void print_prob(FILE * fplog, const char *leg, int n, double *prob);
    void print_transition_matrix(FILE * fplog, int n, int **nmoves, int *nattempt);
    void print_replica_exchange_statistics(FILE * fplog, RecordKeeper * re);
    void print_count(FILE *fplog, const char *leg, int n, int *count);
    void DestroyRecordKeeper();
    bool repl_quantity(vector<Simulation*>* simsRef, RecordKeeper *re);
    void rearrangeByTemperature(int* ind, int* pind, double* temps, int numRepl);
    void rearrangeByTemperature( std::vector<int> & ind, 
                                                        std::vector<int> & pind, 
                                                        double* temps,
                                                        int numRepl);
    void rearrangeByChemPots(int* ind, int* pind, double** chemPots, int numRepl, uint numKinds);
    void rearrangeByChemPots(std::vector<int> & ind, std::vector<int> & pind, double** chemPots, int numRepl, uint numKinds);


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
    RecordKeeper re;
    FILE * fplog;
    Clock * timer;
    bool REMC;

    #if ENSEMBLE == GCMC
    std::map< double, std::vector<int> > temp_map;
    std::map< std::vector<double>, std::vector<int> > mu_map;
    #endif
};

#endif