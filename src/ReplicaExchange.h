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
        multiSimTitle()
    {
    };

    int exchangeInterval; /* Interval in steps at which to attempt exchanges, 0 means no replica exchange */
    int numExchanges;     /* The number of exchanges to attempt at an exchange step */
    int randomSeed;       /* The random seed, -1 means generate a seed */
    std::string multiSimTitle;
};

class ReplicaExchange
{
public:
   // explicit ReplicaExchange(std::vector<Simulation*>*);
    //~ReplicaExchange();

private:

};

#endif