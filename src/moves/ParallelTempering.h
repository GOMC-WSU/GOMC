/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
/*
 * This source code follows the replica exchange design pattern of 
 * GROMACS molecular simulation package, and there is source code within this file
 * that was directly taken and adapted to fit our purposes.
 * citation: GROMACS Abraham, et al. (2015) SoftwareX 1-2 19-25
 */
#ifndef PARALLELTEMPERING_H
#define PARALLELTEMPERING_H

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"
#include <cmath>
#include <mpi.h>
#include "ParallelTemperingMPIMethods.h"

//! Helps cut off probability values.
#define c_probabilityCutoff 100

class ParallelTempering : public MoveBase
{
public:
  ParallelTempering(System &sys, StaticVals const& statV, MultiSim *& multisim, ulong & stepRef);

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();
  //void PrepCFCMC(const uint box);

private:
  void InitRecordKeeper();
  void DestroyRecordKeeper();

  void print_ind(FILE *fplog, const char *leg, int n, int *ind, const bool *bEx);
  void print_prob(FILE *fplog, const char *leg, int n, double *prob);
  double calc_delta(FILE * fplog, bool bPrint, int i, int j, int ip, int jp);
  void print_transition_matrix(FILE * fplog, int n, int **nmoves, int *nattempt);
  void print_count(FILE *fplog, const char *leg, int n, int *count);
  void print_replica_exchange_statistics(FILE * fplog);

  void cyclic_decomposition(const int *destinations, int **cyclic, bool *incycle, const int nrepl, int *nswap);
  void compute_exchange_order(int **cyclic, int **order, const int nrepl, const int maxswap);
  void prepare_to_do_exchange(const int replica_id, int *maxswap, bool *bThisReplicaExchanged);
  void exchange_state(int b);
  void exchange_doubles(int b, double *v, int n);


  // This will require a custom mpi sender
  SystemPotential sysPotNew;
  // These both will need an xyz over mpi sender
  Coordinates newMolsPos;
  COM newCOMs;

  FILE * fplog;

  ulong & step;

  const MoleculeLookup& molLookupRef;

  //! Replica ID
  int       repl;
  //! Total number of replica
  int       nrepl;
  //! Temperature
  double      temp;
  //! Replica exchange type from ere enum
  int       type;
  //! Quantity, e.g. temperature or lambda; first index is ere, second index is replica ID
  double    **q;
  //! Use constant pressure and temperature
  bool  bNPT;
  //! Replica pressures
  double     *pres;
  //! Replica indices
  int      *ind;
  //! Used for keeping track of all the replica swaps
  int      *allswaps;
  //! Replica exchange interval (number of steps)
  const ulong       nst;
  //! Number of exchanges per interval
  int       nex;
  //! Random seed
  int       seed;
  //! Number of even and odd replica change attempts
  int       nattempt[2];
  //! Sum of probabilities
  double     *prob_sum;
  //! Number of moves between replicas i and j
  int     **nmoves;
  //! i-th element of the array is the number of exchanges between replica i-1 and i
  int      *nexchange;

    /*! \brief Helper arrays for replica exchange; allocated here
    * so they don't have to be allocated each time */
  //! \{
  int      *destinations;
  int     **cyclic;
  int     **order;
  int      *tmpswap;
  bool *incycle;
  bool *bEx;

      //! Helper arrays to hold the quantities that are exchanged.
  //! \{
  double  *prob;
  double  *Epot;
  double  *beta;
  double  *Vol;
  double **de;
  //! \}
};

inline ParallelTempering::ParallelTempering(System &sys, StaticVals const &statV, MultiSim *& multisim, ulong & stepRef) :
  MoveBase(sys, statV), molLookupRef(sys.molLookup), step(stepRef), nst(multisim->nst), newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef,statV.mol)
{
  repl = multisim->worldRank;
  nrepl = multisim->worldSize;
  InitRecordKeeper();

  std::string fileName = multisim->pathToReplicaDirectory + "ReplicaLog.dat";
  fplog = fopen(fileName.c_str(), "w");

  fprintf(fplog, "Replica exchange information below: ex and x = exchange, pr = probability\n");

  beta[repl] = statV.forcefield.beta;
  ParallelTemperingMPIMethods::gomc_sumd_comm(nrepl, beta, MPI_COMM_WORLD);
}

inline void ParallelTempering::PrintAcceptKind() {
  printf("%-37s", "% Accepted Parallel Tempering ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::PARALLEL_TEMPERING));
  }
  std::cout << std::endl;
}

inline uint ParallelTempering::Prep(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
  for(int i = 0; i < nrepl; i++){
    Epot[i] = 0.0;
  }
  Epot[repl] = sysPotRef.Total();
  ParallelTemperingMPIMethods::gomc_sumd_comm(nrepl, Epot, MPI_COMM_WORLD);

  return state;
}

inline uint ParallelTempering::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.
  uint state = mv::fail_state::NO_FAIL;
  return state;
}

inline void ParallelTempering::CalcEn()
{
  double delta = 0.0;
  double ediff;
  double probability;
  int     i;
  int     a;
  int     b;
  int     ap;
  int     bp;
  int     m;
  int     tmp;
  bool    bPrint;

  int     *pind     = destinations; /* permuted index */

  // standard nearest neighbor replica exchange /

  /* make a duplicate set of indices for shuffling */
  for (i = 0; i < nrepl; i++)
  {
      pind[i] = ind[i];
  }

  m = (step / nst) % 2;
  for (i = 1; i < nrepl; i++)
  {
      a = ind[i-1];
      b = ind[i];

      bPrint = (repl == a || repl == b);
      if (i % 2 == m)
      {
          delta = calc_delta(fplog, bPrint, a, b, a, b);
          if (delta <= 0)
          {
              // accepted /
              prob[i] = 1;
              bEx[i]  = true;
          }
          else
          {
              if (delta > c_probabilityCutoff)
              {
                  prob[i] = 0;
              }
              else
              {
                  prob[i] = exp(-delta);
              }
              // roll a number to determine if accepted. For now it is superfluous to
              // reset, but just in case we ever add more calls in different branches
              // it is safer to always reset the distribution.
              //uniformRealDist.reset();
              bEx[i] = prng() < prob[i];
          }
          prob_sum[i] += prob[i];

          if (bEx[i])
          {
              // swap these two /
              tmp       = pind[i-1];
              pind[i-1] = pind[i];
              pind[i]   = tmp;
              nexchange[i]++;  /// statistics for back compatibility /
          }
      }
      else
      {
          prob[i] = -1;
          bEx[i]  = false;
      }
  }
  // print some statistics /
  print_ind(fplog, "ex", nrepl, ind, bEx);
  print_prob(fplog, "pr", nrepl, prob);
  fprintf(fplog, "\n");
  nattempt[m]++;
  
  //record which moves were made and accepted 
  for (i = 0; i < nrepl; i++)
  {
      nmoves[ind[i]][pind[i]] += 1;
      nmoves[pind[i]][ind[i]] += 1;
  }  
  fflush(fplog); /* make sure we can see what the last exchange was */

}

inline void ParallelTempering::Accept(const uint rejectState, const uint step)
{
  int j;
  int replica_id = 0;
  int exchange_partner;
  int maxswap = 0;
  /* Number of rounds of exchanges needed to deal with any multiple
    * exchanges. */
  /* Where each replica ends up after the exchange attempt(s). */
  /* The order in which multiple exchanges will occur. */
  bool bThisReplicaExchanged = false;
  replica_id  = repl;
  prepare_to_do_exchange(replica_id, &maxswap, &bThisReplicaExchanged);

  if (bThisReplicaExchanged)
  {
    /* There will be only one swap cycle with standard replica
      * exchange, but there may be multiple swap cycles if we
      * allow multiple swaps. */

    for (j = 0; j < maxswap; j++)
    {
      exchange_partner = order[replica_id][j];

      if (exchange_partner != replica_id)
      {
          fprintf(fplog, "Exchanging %d with %d\n", replica_id, exchange_partner);
          fprintf(fplog, "B4 Exch, My first coords x:%f, y:%f, z:%f\n", coordCurrRef.x[0], coordCurrRef.y[0], coordCurrRef.z[0]);
          fflush(fplog);
          exchange_state(exchange_partner);
          fprintf(fplog, "AF Exch, My first coords x:%f, y:%f, z:%f\n", coordCurrRef.x[0], coordCurrRef.y[0], coordCurrRef.z[0]);
          fflush(fplog);
      }
    }
  }

  bool result = bThisReplicaExchanged;

  if(result){
    cellList.GridAll(boxDimRef, coordCurrRef, molLookupRef);

  }

  uint mkTot = molLookupRef.GetNumCanMoveKind();

  for(uint b = 0; b < BOX_TOTAL; b++) {
    for (uint mk = 0; mk < mkTot; ++mk){
      moveSetRef.Update(mv::PARALLEL_TEMPERING, result, step, b, mk);
    }
  }
}

void ParallelTempering::InitRecordKeeper()
{
  nattempt[0] = 0;
  nattempt[1] = 0;
  prob_sum = new double[nrepl]();
  nexchange = new int[nrepl]();
  nmoves = new int*[nrepl];
  for (int i = 0; i < nrepl; i++) {
    nmoves[i] = new int[nrepl]();
  }
  ind = new int[nrepl]();

  for (int i = 0; i < nrepl; i++)
  {
      ind[i] = i;
  }

  #if ENSEMBLE == NPT
    bNPT = true;
  #else
    bNPT = false;
  #endif

  destinations = new int[nrepl];
  incycle = new bool[nrepl];
  tmpswap = new int[nrepl];
  cyclic = new int*[nrepl];
  order = new int*[nrepl];
  for (int i = 0; i < nrepl; i++)
  {
      cyclic[i] = new int[nrepl+1];
      order[i] = new int[nrepl];
  }
  /* allocate space for the functions storing the data for the replicas */
  /* not all of these arrays needed in all cases, but they don't take
      up much space, since the max size is nrepl**2 */
  prob = new double[nrepl]();
  bEx =  new bool[nrepl]();
  beta = new double[nrepl]();
  Vol = new double[nrepl]();
  Epot = new double[nrepl]();
  de = new double*[nrepl];
  for (int i = 0; i < nrepl; i++)
  {
    de[i] = new double[nrepl];
  }
}

void ParallelTempering::DestroyRecordKeeper(){
  for (int i = 0; i < nrepl; i++)
  {
    delete[] de[i];
  }
  delete[] de;
  delete[] Epot;
  delete[] Vol;
  delete[] beta;
  delete[] bEx;
  delete[] prob;
  for (int i = 0; i < nrepl; i++) 
  {
    delete[] order[i];
  }
  for (int i = 0; i < nrepl+1; i++) {
    delete[] cyclic[i];
  }
  delete[] order;
  delete[] cyclic;
  delete[] tmpswap;
  delete[] incycle;
  delete[] destinations;
  delete[] ind;
  for (int i = 0; i < nrepl; i++) {
    delete[] nmoves[i];
  }
  delete[] nmoves;
  delete[] nexchange;
  delete[] prob_sum;
}

void ParallelTempering::print_ind(FILE *fplog, const char *leg, int n, int *ind, const bool *bEx)
{
    int i;

    fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %c %2d", (bEx != NULL && bEx[i]) ? 'x' : ' ', ind[i]);
    }
    fprintf(fplog, "\n");
}

void ParallelTempering::print_prob(FILE *fplog, const char *leg, int n, double *prob)
{
    int  i;
    char buf[8];

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf, "%4.2f", prob[i]);
            fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog, "     ");
        }
    }
    fprintf(fplog, "\n");
}

double ParallelTempering::calc_delta(FILE * fplog, bool bPrint, int i, int j, int ip, int jp){

  double delta;

  #if ENSEMBLE == NPT || ENSEMBLE == NVT || ENSEMBLE == GCMC

  double U_i = Epot[i];
  double U_j = Epot[j];
  double beta_i = beta[i];
  double beta_j = beta[j];

  double ediff = U_j - U_i;
  delta = -(beta_j - beta_i)*ediff;

  #endif

  // GROMACS Abraham, et al. (2015) SoftwareX 1-2 19-25 /
  #if ENSEMBLE == NPT || ENSEMBLE == NVT
  if(bPrint)
    fprintf(fplog, "Repl %d <-> %d  dE = %10.3e, dBeta = %10.3e, dBeta_dE = %10.3e \n", i, j, ediff, beta_j - beta_i, delta);
  #endif
/*
  //  GROMACS Abraham, et al. (2015) SoftwareX 1-2 19-25 /
  #if ENSEMBLE == NPT
    double pres_i = (*simsRef)[ip]->getPressure();
    double vol_i = (*simsRef)[i]->getVolume();
    double pres_j = (*simsRef)[jp]->getPressure();
    double vol_j = (*simsRef)[j]->getVolume();
    double dpV = (beta_i * pres_i - beta_j * pres_j) * (vol_j - vol_i);
    fprintf(fplog, "  dpV = %10.3e  dBeta_dE + dpV = %10.3e\n", dpV, delta + dpV);
    delta += dpV;
  #endif

  //  Faller, Roland, Qiliang Yan, and Juan J. de Pablo. "Multicanonical parallel tempering." The Journal of chemical physics 116.13 (2002): 5419-5423.
   //   Eq (3)/
  #if ENSEMBLE == GCMC
  // TODO
  // / Add support for fugacity
   //
    double deltaBetaMuN   = 0;

    for (uint index = 0; index < ((*simsRef)[i]->getSystem())->molLookup.GetNumKind(); index++){
      deltaBetaMuN -= ((beta_i * (*simsRef)[ip]->getChemicalPotential(index) -
        beta_j * (*simsRef)[jp]->getChemicalPotential(index) ) * (
         (*simsRef)[j]->getNumOfParticles(index) - (*simsRef)[i]->getNumOfParticles(index)));
    }
    fprintf(fplog, "Repl %d <-> %d  dE = %10.3e, dBeta = %10.3e, dBeta_dE = %10.3e \n", i, j, ediff, beta_j - beta_i, delta);
    fprintf(fplog, "  dMuN = %10.3e  dBeta_dE + dMuN = %10.3e\n", deltaBetaMuN, delta + deltaBetaMuN);
    delta += deltaBetaMuN;
  #endif

  //  Ortiz et al Chemical physics letters 368.3-4 (2003): 452-457.
  //    Eq (8) /
  #if ENSEMBLE == GEMC

  delta = 0;
  double dpV = 0;
  double prob = 1.0;
  double individual_delta = 0.0;
  for (uint box_index = 0; box_index < BOX_TOTAL; ++box_index) { 
      delta += -(((*simsRef)[jp]->getBeta() - (*simsRef)[ip]->getBeta())*
        ((*simsRef)[j]->getEpotBox(box_index) - (*simsRef)[i]->getEpotBox(box_index)));
      individual_delta =  -(((*simsRef)[jp]->getBeta() - (*simsRef)[ip]->getBeta())*
        ((*simsRef)[j]->getEpotBox(box_index) - (*simsRef)[i]->getEpotBox(box_index)));
      fprintf(fplog, "Repl %d <-> %d  BOX_%u : probability = %10.3e \n", i, j, box_index, min(1.0, exp(-individual_delta)));
      prob *= min(1.0, exp(-individual_delta));
      if ((*simsRef)[i]->getKindOfGEMC() == mv::GEMC_NPT) {
        dpV += ((*simsRef)[i]->getBeta() * (*simsRef)[ip]->getPressure() - 
          (*simsRef)[j]->getBeta() * (*simsRef)[jp]->getPressure()) * 
          ((*simsRef)[j]->getVolume(box_index) - (*simsRef)[i]->getVolume(box_index));
      }
  }
  fprintf(fplog, "Repl %d <-> %d  BOX_TOTAL : probability by a0*a1 = %10.3e \n", i, j, prob);
  fprintf(fplog, "Repl %d <-> %d  BOX_TOTAL : probability by prob sum of deltas = %10.3e \n", i, j, min(1.0, exp(-delta)));

  fprintf(fplog, "Repl %d <-> %d  BOX_TOTAL : dBeta_dE = %10.3e \n", i, j, delta);
  if ((*simsRef)[i]->getKindOfGEMC() == mv::GEMC_NPT) {
    fprintf(fplog, "BOX_TOTAL  dpV = %10.3e  dBeta_dE + dpV = %10.3e\n", dpV, delta + dpV);
  }
  delta += dpV; 

  #endif
*/
  return delta;
}

void ParallelTempering::print_transition_matrix(FILE * fplog, int n, int **nmoves, int *nattempt)
{
     int   i, j, ntot;
    float Tprint;

    ntot = nattempt[0] + nattempt[1];
    fprintf(fplog, "\n");
    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "    ");  /* put the title closer to the center */
    }
    fprintf(fplog, "Empirical Transition Matrix\n");

    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%8d", (i+1));
    }
    fprintf(fplog, "\n");

    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "Repl");
        for (j = 0; j < n; j++)
        {
            Tprint = 0.0;
            if (nmoves[i][j] > 0)
            {
                Tprint = nmoves[i][j]/(2.0*ntot);
            }
            fprintf(fplog, "%8.4f", Tprint);
        }
        fprintf(fplog, "%3d\n", i);
    }
}

void ParallelTempering::print_replica_exchange_statistics(FILE * fplog){
    int  i;

    fprintf(fplog, "\nReplica exchange statistics\n");

   // if (nex == 0)
   // {
        fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
                nattempt[0]+nattempt[1], nattempt[1], nattempt[0]);

        fprintf(fplog, "Repl  average probabilities:\n");
        for (i = 1; i < nrepl; i++)
        {
            if (nattempt[i%2] == 0)
            {
                prob[i] = 0;
            }
            else
            {
                prob[i] =  prob_sum[i]/nattempt[i%2];
            }
        }
        print_ind(fplog, "", nrepl, ind, NULL);
        print_prob(fplog, "", nrepl, prob);

        fprintf(fplog, "Repl  number of exchanges:\n");
        print_ind(fplog, "", nrepl, ind, NULL);
        print_count(fplog, "", nrepl, nexchange);

        fprintf(fplog, "Repl  average number of exchanges:\n");
        for (i = 1; i < nrepl; i++)
        {
            if (nattempt[i%2] == 0)
            {
                prob[i] = 0;
            }
            else
            {
                prob[i] =  ((double)nexchange[i])/nattempt[i%2];
            }
        }
        print_ind(fplog, "", nrepl, ind, NULL);
        print_prob(fplog, "", nrepl, prob);

        fprintf(fplog, "\n");
    //}
    /* print the transition matrix */
    print_transition_matrix(fplog, nrepl, nmoves, nattempt);
}

void ParallelTempering::print_count(FILE *fplog, const char *leg, int n, int *count){
    int i;

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %4d", count[i]);
    }
    fprintf(fplog, "\n");
}

void
ParallelTempering::cyclic_decomposition(const int *destinations,
                     int      **cyclic,
                     bool  *incycle,
                     const int  nrepl,
                     int       *nswap)
{

    int i, j, c, p;
    int maxlen = 1;
    for (i = 0; i < nrepl; i++)
    {
        incycle[i] = false;
    }
    for (i = 0; i < nrepl; i++)  /* one cycle for each replica */
    {
        if (incycle[i])
        {
            cyclic[i][0] = -1;
            continue;
        }
        cyclic[i][0] = i;
        incycle[i]   = true;
        c            = 1;
        p            = i;
        for (j = 0; j < nrepl; j++) /* Epotly all cycles are part, but we will break first */
        {
            p = destinations[p];    /* start permuting */
            if (p == i)
            {
                cyclic[i][c] = -1;
                if (c > maxlen)
                {
                    maxlen = c;
                }
                break; /* we've reached the original element, the cycle is complete, and we marked the end. */
            }
            else
            {
                cyclic[i][c] = p;  /* each permutation gives a new member of the cycle */
                incycle[p]   = true;
                c++;
            }
        }
    }
    *nswap = maxlen - 1;

#ifndef NDEBUG

    for (i = 0; i < nrepl; i++)
    {
        fprintf(fplog, "Cycle %d:", i);
        for (j = 0; j < nrepl; j++)
        {
            if (cyclic[i][j] < 0)
            {
                break;
            }
            fprintf(fplog, "%2d", cyclic[i][j]);
        }
        fprintf(fplog, "\n");
    }
    fflush(fplog);
    
#endif
}

void
ParallelTempering::compute_exchange_order(int     **cyclic,
                       int     **order,
                       const int nrepl,
                       const int maxswap)
{
    int i, j;

    for (j = 0; j < maxswap; j++)
    {
        for (i = 0; i < nrepl; i++)
        {
            if (cyclic[i][j+1] >= 0)
            {
                order[cyclic[i][j+1]][j] = cyclic[i][j];
                order[cyclic[i][j]][j]   = cyclic[i][j+1];
            }
        }
        for (i = 0; i < nrepl; i++)
        {
            if (order[i][j] < 0)
            {
                order[i][j] = i; /* if it's not exchanging, it should stay this round*/
            }
        }
    }

#ifndef NDEBUG

        fprintf(fplog, "Replica Exchange Order\n");
        for (i = 0; i < nrepl; i++)
        {
            fprintf(fplog, "Replica %d:", i);
            for (j = 0; j < maxswap; j++)
            {
                if (order[i][j] < 0)
                {
                    break;
                }
                fprintf(fplog, "%2d", order[i][j]);
            }
            fprintf(fplog, "\n");
        }
        fflush(fplog);
#endif
}

void
ParallelTempering::prepare_to_do_exchange(const int           replica_id,
                                                int                *maxswap,
                                                bool           *bThisReplicaExchanged)
{
    int i, j;
    /* Hold the cyclic decomposition of the (multiple) replica
     * exchange. */
    bool bAnyReplicaExchanged = false;
    *bThisReplicaExchanged = false;

    for (i = 0; i < nrepl; i++)
    {
        if (destinations[i] != ind[i])
        {
            /* only mark as exchanged if the index has been shuffled */
            bAnyReplicaExchanged = true;
            break;
        }
    }
    if (bAnyReplicaExchanged)
    {
        /* reinitialize the placeholder arrays */
        for (i = 0; i < nrepl; i++)
        {
            for (j = 0; j < nrepl; j++)
            {
                cyclic[i][j] = -1;
                order[i][j]  = -1;
            }
        }

        /* Identify the cyclic decomposition of the permutation (very
         * fast if neighbor replica exchange). */
        cyclic_decomposition(destinations, cyclic, incycle, nrepl, maxswap);

        /* Now translate the decomposition into a replica exchange
         * order at each step. */
        compute_exchange_order(cyclic, order, nrepl, *maxswap);

        /* Did this replica do any exchange at any point? */
        for (j = 0; j < *maxswap; j++)
        {
            if (replica_id != order[replica_id][j])
            {
                *bThisReplicaExchanged = true;
                break;
            }
        }
    }
}

void ParallelTempering::exchange_state(int b)
{

  //sysPotRef = sysPotNew;
  //swap(coordCurrRef, newMolsPos);
  //swap(comCurrRef, newCOMs);
  //update reciprocate value
  //calcEwald->UpdateRecip(bPick);
  exchange_doubles(b, &sysPotRef.totalEnergy.correction, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.inter, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.intraBond, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.intraNonbond, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.real, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.recip, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.self, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.tc, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.total, 1);
  exchange_doubles(b, &sysPotRef.totalEnergy.totalElect, 1);
  
  exchange_doubles(b, coordCurrRef.x, coordCurrRef.Count());
  exchange_doubles(b, coordCurrRef.y, coordCurrRef.Count());
  exchange_doubles(b, coordCurrRef.z, coordCurrRef.Count());
  
  exchange_doubles(b, comCurrRef.x, comCurrRef.Count());
  exchange_doubles(b, comCurrRef.y, comCurrRef.Count());
  exchange_doubles(b, comCurrRef.z, comCurrRef.Count());
}

void ParallelTempering::exchange_doubles(int b, double *v, int n)
{
    double *buf;
    int   i;

    if (v)
    {
      buf =  new double[n];
        /*
           MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(double), MPI_BYTE, b, 0,
                      MPI_COMM_WORLD, &mpi_req);
            MPI_Recv(buf, n*sizeof(double), MPI_BYTE, b, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
        for (i = 0; i < n; i++)
        {
            //fprintf(fplog, "%f\n", buf[i]);
            //fflush(fplog);
            v[i] = buf[i];
        }
        delete[] buf;
    }
}
#endif
