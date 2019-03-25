/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef REPLICAEXCHANGE_H
#define REPLICAEXCHANGE_H

//Member vars
/* Courtest of GROMACS */
/* The parameters for the replica exchange algorithm */
struct ReplicaExchangeParameters
{
    ReplicaExchangeParameters() :
        exchangeInterval(0),
        numExchanges(0),
        randomSeed(-1),
        multiSimTitle("Replica_Exchange_Simulation")
    {
    };

    ulong exchangeInterval; /* Interval in steps at which to attempt exchanges, 0 means no replica exchange */
    int numExchanges;     /* The number of exchanges to attempt at an exchange step */
    int randomSeed;       /* The random seed, -1 means generate a seed */
    std::string multiSimTitle;
};
struct RecordKeeper
{
    int         nrepl;
    int       nattempt[2]; /* number of even and odd replica change attempts */
    bool        *bEx;   
    double   *prob;    /* probabilities */
    double   *prob_sum;    /* sum of probabilities */
    int     **nmoves;      /* number of moves between replicas i and j */
    int      *nexchange;   /* i-th element of the array is the number of exchanges between replica i-1 and i */
    int     *ind;       /* An array of indices illustrating where our temps 0..N are orderwise */
    int     *pind;       /* A Permuttable Map of indices illustrating where our temps 0..N are orderwise */
    map<int, int> temperatureToIndex;
};

#endif