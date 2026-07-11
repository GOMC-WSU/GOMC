/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011-2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief Implements the replica exchange routines.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */

/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in License.txt
along with this program. also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingUtilities.h"

#if GOMC_LIB_MPI

ParallelTemperingUtilities::ParallelTemperingUtilities(
    MultiSim const *const &multisim, System &sys, StaticVals const &statV,
    ulong parallelTempFreq, ulong parallelTemperingAttemptsPerExchange)
    : ms(multisim), fplog(multisim->fplog), sysPotRef(sys.potential),
      parallelTempFreq(parallelTempFreq),
      parallelTemperingAttemptsPerExchange(
          parallelTemperingAttemptsPerExchange),
      prng(*sys.prngParallelTemp),
      newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
      newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol) {
  if (parallelTemperingAttemptsPerExchange > 1) {
    bMultiEx = true;
    fprintf(fplog,
            "Using multiple random switch exchange. Number of Attempts per "
            "Exchange : %lu\n",
            parallelTemperingAttemptsPerExchange);
  } else {
    bMultiEx = false;
    fprintf(fplog, "Using standard nearest neighbor exchange.\n");
  }

#if BOX_TOTAL == 1
  global_energies.resize(ms->worldSize, 0.0);
#else
  global_energies.resize(2, std::vector<double>(ms->worldSize, 0.0));
#endif
  global_betas.resize(ms->worldSize, 0.0);
  ind.resize(ms->worldSize, 0);
  pind.resize(ms->worldSize, 0);
  exchangeProbabilities.resize(ms->worldSize, 0.0);
  prob_sum.resize(ms->worldSize, 0.0);
  exchangeResults.resize(ms->worldSize, false);
  nexchange.resize(ms->worldSize, 0);
  allswaps.resize(ms->worldSize, 0);
  tempswap.resize(ms->worldSize, 0);
  nattempt.resize(2, 0);
  nmoves.resize(ms->worldSize, std::vector<int>(ms->worldSize, 0));
  global_betas[ms->worldRank] = statV.forcefield.beta;
  cyclic.resize(ms->worldSize, std::vector<int>(ms->worldSize + 1, -1));
  order.resize(ms->worldSize, std::vector<int>(ms->worldSize, -1));
  incycle.resize(ms->worldSize, false);

  MPI_Allreduce(MPI_IN_PLACE, &global_betas[0], ms->worldSize, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  for (int i = 0; i < ms->worldSize; i++) {
    ind[i] = i;
    allswaps[i] = i;
  }
}

void ParallelTemperingUtilities::evaluateExchangeCriteria(ulong step) {
  // for (int i = 0; i < 2; i++){
  //  std::cout << "Before fill : energy[" << i << "] : " << global_energies[i]
  //  << std::endl;
  //}
#if BOX_TOTAL == 1
  std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
#else
  std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
  std::memset(&global_energies[1], 0, ms->worldSize * sizeof(double));
#endif

  for (int i = 0; i < ms->worldSize; i++) {
    pind[i] = ind[i];
  }

  // for (int i = 0; i < 2; i++){
  //  std::cout << "After fill : energy[" << i << "] : " << global_energies[i]
  //  << std::endl;
  //}

#if BOX_TOTAL == 1
  //        global_energies[ms->worldRank] = sysPotRef.boxEnergy[0].total;
  global_energies[ms->worldRank] = sysPotRef.totalEnergy.total;

  // for (int i = 0; i < 2; i++){
  //  std::cout << "After set local energy : energy[" << i << "] : " <<
  //  global_energies[i] << std::endl;
  //}

  MPI_Allreduce(MPI_IN_PLACE, &global_energies[0], ms->worldSize, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  // for (int i = 0; i < 2; i++){
  //   std::cout << "After allreduce : energy[" << i << "] : " <<
  //   global_energies[i] << std::endl;
  // }

#else

  global_energies[0][ms->worldRank] = sysPotRef.boxEnergy[0].total;
  global_energies[1][ms->worldRank] = sysPotRef.boxEnergy[1].total;

  MPI_Allreduce(MPI_IN_PLACE, &global_energies[0], ms->worldSize * 2,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // MPI_Allreduce(MPI_IN_PLACE, &energies[1], ms->worldSize, MPI_DOUBLE,
  // MPI_SUM,
  //    MPI_COMM_WORLD);

#endif
  double uBoltz;
  bool bPrint = false;
  double printRecord;
  if (bMultiEx) {
    int i0, i1, a, b, ap, bp;

    /* multiple random switch exchange */
    int nself = 0;

    for (int i = 0; i < parallelTemperingAttemptsPerExchange + nself; i++) {
      i0 = prng.randIntExc(ms->worldSize);
      i1 = prng.randIntExc(ms->worldSize);
      if (i0 == i1) {
        nself++;
        continue; /* self-exchange, back up and do it again */
      }

      a = ind[i0]; /* what are the indices of these states? */
      b = ind[i1];
      ap = pind[i0];
      bp = pind[i1];

      bPrint = false; /* too noisy */
                      /* calculate the energy difference */
#if ENSEMBLE == NVT
      uBoltz = exp((global_betas[b] - global_betas[a]) *
                   (global_energies[bp] - global_energies[ap]));
#endif
      /* we actually only use the first space in the prob and bEx array,
         since there are actually many switches between pairs. */
      exchangeProbabilities[0] = std::min(uBoltz, 1.0);
      exchangeResults[0] = (printRecord = prng()) < uBoltz;
      // std::cout << "Swapping repl " << i-1 << " and repl " << i << " uBoltz
      // :" << uBoltz << "prng : " << printRecord << std::endl;
      prob_sum[0] += exchangeProbabilities[0];
      if (exchangeResults[0]) {
        /* swap these two */
        int tmp = pind[i0];
        pind[i0] = pind[i1];
        pind[i1] = tmp;
      }
    }
    nattempt[0]++; /* keep track of total permutation trials here */
    print_allswitchind(fplog, ms->worldSize, pind, allswaps, tempswap);

  } else {
    int parity = step / parallelTempFreq % 2;

    for (int i = 1; i < ms->worldSize; i++) {
      bPrint = ms->worldRank == i || ms->worldRank == i - 1;
      if (i % 2 == parity) {
#if ENSEMBLE == NVT
        uBoltz = exp((global_betas[i] - global_betas[i - 1]) *
                     (global_energies[i] - global_energies[i - 1]));
#endif
        exchangeProbabilities[i] = std::min(uBoltz, 1.0);
        exchangeResults[i] = (printRecord = prng()) < uBoltz;
        // std::cout << "Swapping repl " << i-1 << " and repl " << i << " uBoltz
        // :" << uBoltz << "prng : " << printRecord << std::endl;
        prob_sum[i] += exchangeProbabilities[i];
        if (exchangeResults[i]) {
          /* swap these two */
          int tmp = pind[i - 1];
          pind[i - 1] = pind[i];
          pind[i] = tmp;
          nexchange[i]++;
        }
      } else {
        exchangeResults[i] = false;
        exchangeProbabilities[i] = 0.0;
      }
    }
    print_ind(fplog, "ex", ms->worldSize, ind, exchangeResults);
    print_prob(fplog, "pr", ms->worldSize, exchangeProbabilities);
    fprintf(fplog, "\n");
    nattempt[parity]++;
  }

  for (int i = 0; i < ms->worldSize; i++) {
    nmoves[ind[i]][pind[i]] += 1;
    nmoves[pind[i]][ind[i]] += 1;
  }
}

void ParallelTemperingUtilities::prepareToDoExchange(
    const int replica_id, int *maxswap, bool *bThisReplicaExchanged) {
  int i, j;
  /* Hold the cyclic decomposition of the (multiple) replica
   * exchange. */
  bool bAnyReplicaExchanged = false;
  *bThisReplicaExchanged = false;

  for (i = 0; i < ms->worldSize; i++) {
    if (pind[i] != ind[i]) {
      /* only mark as exchanged if the index has been shuffled */
      bAnyReplicaExchanged = true;
      break;
    }
  }
  if (bAnyReplicaExchanged) {
    /* reinitialize the placeholder arrays */
    for (i = 0; i < ms->worldSize; i++) {
      for (j = 0; j < ms->worldSize; j++) {
        cyclic[i][j] = -1;
        order[i][j] = -1;
      }
    }

    /* Identify the cyclic decomposition of the permutation (very
     * fast if neighbor replica exchange). */
    cyclicDecomposition(pind, cyclic, incycle, ms->worldSize, maxswap);

    /* Now translate the decomposition into a replica exchange
     * order at each step. */
    computeExchangeOrder(cyclic, order, ms->worldSize, *maxswap);

    /* Did this replica do any exchange at any point? */
    for (j = 0; j < *maxswap; j++) {
      if (replica_id != order[replica_id][j]) {
        *bThisReplicaExchanged = true;
        break;
      }
    }
  }
}

void ParallelTemperingUtilities::cyclicDecomposition(
    const std::vector<int> &destinations, std::vector<std::vector<int>> &cyclic,
    std::vector<bool> &incycle, const int nrepl, int *nswap) {
  int i, j, c, p;
  int maxlen = 1;
  for (i = 0; i < nrepl; i++) {
    incycle[i] = false;
  }
  for (i = 0; i < nrepl; i++) { /* one cycle for each replica */
    if (incycle[i]) {
      cyclic[i][0] = -1;
      continue;
    }
    cyclic[i][0] = i;
    incycle[i] = true;
    c = 1;
    p = i;
    for (j = 0; j < nrepl;
         j++) { /* potentially all cycles are part, but we will break first */
      p = destinations[p]; /* start permuting */
      if (p == i) {
        cyclic[i][c] = -1;
        if (c > maxlen) {
          maxlen = c;
        }
        break; /* we've reached the original element, the cycle is complete, and
                  we marked the end. */
      } else {
        cyclic[i][c] = p; /* each permutation gives a new member of the cycle */
        incycle[p] = true;
        c++;
      }
    }
  }
  *nswap = maxlen - 1;

#ifndef NDEBUG
  for (i = 0; i < nrepl; i++) {
    fprintf(fplog, "Cycle %d:", i);
    for (j = 0; j < nrepl; j++) {
      if (cyclic[i][j] < 0) {
        break;
      }
      fprintf(fplog, "%2d", cyclic[i][j]);
    }
    fprintf(fplog, "\n");
  }
  fflush(fplog);
#endif
}

void ParallelTemperingUtilities::computeExchangeOrder(
    std::vector<std::vector<int>> &cyclic, std::vector<std::vector<int>> &order,
    const int nrepl, const int maxswap) {
  int i, j;
  for (j = 0; j < maxswap; j++) {
    for (i = 0; i < nrepl; i++) {
      if (cyclic[i][j + 1] >= 0) {
        order[cyclic[i][j + 1]][j] = cyclic[i][j];
        order[cyclic[i][j]][j] = cyclic[i][j + 1];
      }
    }
    for (i = 0; i < nrepl; i++) {
      if (order[i][j] < 0) {
        order[i][j] = i; /* if it's not exchanging, it should stay this round*/
        fprintf(fplog, "order[%d][%d] = %d\n", i, j, order[i][j]);

        fflush(fplog);
      }
    }
  }

#ifndef NDEBUG

  fprintf(fplog, "Replica Exchange Order\n");
  for (i = 0; i < nrepl; i++) {
    fprintf(fplog, "Replica %d:", i);
    for (j = 0; j < maxswap; j++) {
      if (order[i][j] < 0) {
        break;
      }
      fprintf(fplog, "%2d", order[i][j]);
    }
    fprintf(fplog, "\n");
  }
  fflush(fplog);
#endif
}

void ParallelTemperingUtilities::conductExchanges(
    Coordinates &coordCurrRef, COM &comCurrRef, MultiSim const *const &ms,
    const int &maxSwap, const bool &bThisReplicaExchanged) {
  int exchangePartner;
  int replicaID = ms->worldRank;

  if (bThisReplicaExchanged) {
    for (int j = 0; j < maxSwap; j++) {
      exchangePartner = order[replicaID][j];

      if (exchangePartner != replicaID) {
        /* Exchange the global states between the master nodes */
#ifndef NDEBUG
        fprintf(fplog, "Exchanging %d with %d\n", replicaID, exchangePartner);
#endif

        newMolsPos = coordCurrRef;
        newCOMs = comCurrRef;

        exchangePositionsNonBlocking(newMolsPos, ms, exchangePartner);
        exchangeCOMsNonBlocking(newCOMs, ms, exchangePartner);

        swap(coordCurrRef, newMolsPos);
        swap(comCurrRef, newCOMs);
      }
    }
  }
}

void ParallelTemperingUtilities::exchangePositionsNonBlocking(
    Coordinates &myPos, MultiSim const *const &multisim, int exchangePartner) {
  XYZArray buffer(myPos);

  MPI_Request mpi_req;

  MPI_Isend(buffer.x, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos.x, myPos.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.y, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos.y, myPos.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.z, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos.z, myPos.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
}

void ParallelTemperingUtilities::exchangeCOMsNonBlocking(
    COM &myCOMs, MultiSim const *const &multisim, int exchangePartner) {
  XYZArray buffer(myCOMs);

  MPI_Request mpi_req;

  MPI_Isend(buffer.x, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs.x, myCOMs.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.y, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs.y, myCOMs.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.z, buffer.Count() * sizeof(double), MPI_BYTE,
            exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs.z, myCOMs.Count() * sizeof(double), MPI_BYTE, exchangePartner,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
}

void ParallelTemperingUtilities::exchangePositions(
    Coordinates &myPos, MultiSim const *const &multisim, int exchangePartner,
    bool leader) {
  XYZArray buffer(myPos);

  // if im 0, im the follower and i get 1 as a
  if (leader) {
    MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
  if (leader) {
    MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
  if (leader) {
    MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
}

void ParallelTemperingUtilities::exchangeCOMs(COM &myCOMs,
                                              MultiSim const *const &multisim,
                                              int exchangePartner,
                                              bool leader) {
  XYZArray buffer(myCOMs);

  if (leader) {
    MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
  if (leader) {
    MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
  if (leader) {
    MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }
}

void ParallelTemperingUtilities::exchangeCellLists(
    CellList &myCellList, MultiSim const *const &multisim, int exchangePartner,
    bool leader) {
  CellList buffer(myCellList);
  std::cout << " buffer list size : " << buffer.list.size()
            << " myCellList list Size : " << myCellList.list.size()
            << std::endl;
  std::cout << " buffer head size : " << buffer.head[0].size()
            << " myCellList head Size : " << myCellList.head[0].size()
            << std::endl;

  // if im 0, im the follower and i get 1 as a

  if (leader) {
    MPI_Send(&buffer.list[0], buffer.list.size(), MPI_INT, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&myCellList.list[0], myCellList.list.size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&myCellList.list[0], myCellList.list.size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.list[0], buffer.list.size(), MPI_INT, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.head[0][0], buffer.head[0].size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD);
    MPI_Recv(&myCellList.head[0][0], myCellList.head[0].size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&myCellList.head[0][0], myCellList.head[0].size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.head[0][0], buffer.head[0].size(), MPI_INT,
             exchangePartner, 0, MPI_COMM_WORLD);
  }
}

void ParallelTemperingUtilities::exchangePotentials(
    SystemPotential &mySystemPotential, MultiSim const *const &multisim,
    int exchangePartner, bool leader) {
  SystemPotential buffer(mySystemPotential);

  // if im 0, im the follower and i get 1 as a

  if (leader) {
    MPI_Send(&buffer.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.total, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.total, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.totalElect, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.totalElect, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.correction, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.correction, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.inter, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.inter, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.intraBond, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.intraBond, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.intraNonbond, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.intraNonbond, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.real, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.real, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.recip, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.recip, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.self, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.self, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalEnergy.tailCorrection, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalEnergy.tailCorrection, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalEnergy.tailCorrection, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalEnergy.tailCorrection, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
  }
  for (int b = 0; b < BOX_TOTAL; b++) {
    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].total, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].total, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].totalElect, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].totalElect, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].correction, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].correction, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].inter, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].inter, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].intraBond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].intraBond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].real, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].real, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].recip, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].recip, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].self, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].self, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxEnergy[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxEnergy[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxEnergy[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.totalEnergy.tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
    }
  }
}

void ParallelTemperingUtilities::exchangeVirials(
    SystemPotential &mySystemPotential, MultiSim const *const &multisim,
    int exchangePartner, bool leader) {
  SystemPotential buffer(mySystemPotential);

  // if im 0, im the follower and i get 1 as a

  if (leader) {
    MPI_Send(&buffer.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.total, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.total, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.totalElect, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.totalElect, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.correction, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.correction, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.inter, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.inter, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.real, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.real, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.recip, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.recip, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.self, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.self, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
             MPI_COMM_WORLD);
  }

  if (leader) {
    MPI_Send(&buffer.totalVirial.tailCorrection, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
    MPI_Recv(&mySystemPotential.totalVirial.tailCorrection, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&mySystemPotential.totalVirial.tailCorrection, 1, MPI_DOUBLE,
             exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&buffer.totalVirial.tailCorrection, 1, MPI_DOUBLE, exchangePartner,
             0, MPI_COMM_WORLD);
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (leader) {
        MPI_Send(&buffer.totalVirial.totalTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.totalTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&mySystemPotential.totalVirial.totalTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.totalTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
      }

      if (leader) {
        MPI_Send(&buffer.totalVirial.interTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.interTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&mySystemPotential.totalVirial.interTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.interTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
      }

      if (leader) {
        MPI_Send(&buffer.totalVirial.realTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.realTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&mySystemPotential.totalVirial.realTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.realTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
      }

      if (leader) {
        MPI_Send(&buffer.totalVirial.recipTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.recipTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&mySystemPotential.totalVirial.recipTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.recipTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
      }

      if (leader) {
        MPI_Send(&buffer.totalVirial.corrTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.corrTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      } else {
        MPI_Recv(&mySystemPotential.totalVirial.corrTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.corrTens[i][j], 1, MPI_DOUBLE,
                 exchangePartner, 0, MPI_COMM_WORLD);
      }
    }
  }

  for (int b = 0; b < BOX_TOTAL; b++) {
    if (leader) {
      MPI_Send(&buffer.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].total, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].total, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].totalElect, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].totalElect, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].correction, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].correction, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner,
               0, MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].inter, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].inter, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].real, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].real, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].recip, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].recip, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].self, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].self, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
               MPI_COMM_WORLD);
    }

    if (leader) {
      MPI_Send(&buffer.boxVirial[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
      MPI_Recv(&mySystemPotential.boxVirial[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(&mySystemPotential.boxVirial[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&buffer.boxVirial[b].tailCorrection, 1, MPI_DOUBLE,
               exchangePartner, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (leader) {
          MPI_Send(&buffer.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
          MPI_Recv(&mySystemPotential.boxVirial[b].totalTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(&mySystemPotential.boxVirial[b].totalTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Send(&buffer.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
        }

        if (leader) {
          MPI_Send(&buffer.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
          MPI_Recv(&mySystemPotential.boxVirial[b].interTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(&mySystemPotential.boxVirial[b].interTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Send(&buffer.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
        }

        if (leader) {
          MPI_Send(&buffer.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
          MPI_Recv(&mySystemPotential.boxVirial[b].realTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(&mySystemPotential.boxVirial[b].realTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Send(&buffer.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
        }

        if (leader) {
          MPI_Send(&buffer.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
          MPI_Recv(&mySystemPotential.boxVirial[b].recipTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(&mySystemPotential.boxVirial[b].recipTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Send(&buffer.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
        }

        if (leader) {
          MPI_Send(&buffer.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
          MPI_Recv(&mySystemPotential.boxVirial[b].corrTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        } else {
          MPI_Recv(&mySystemPotential.boxVirial[b].corrTens[i][j], 1,
                   MPI_DOUBLE, exchangePartner, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          MPI_Send(&buffer.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE,
                   exchangePartner, 0, MPI_COMM_WORLD);
        }
      }
    }
  }
}

void ParallelTemperingUtilities::print_ind(FILE *fplog, const char *leg, int n,
                                           const std::vector<int> &ind,
                                           const std::vector<bool> &bEx) {
  int i;

  fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
  for (i = 1; i < n; i++) {
    fprintf(fplog, " %c %2d", (bEx.empty() != true && bEx[i]) ? 'x' : ' ',
            ind[i]);
  }
  fprintf(fplog, "\n");
}

void ParallelTemperingUtilities::print_prob(FILE *fplog, const char *leg, int n,
                                            const std::vector<double> &prob) {
  int i;
  char buf[8];

  fprintf(fplog, "Repl %2s ", leg);
  for (i = 1; i < n; i++) {
    if (prob[i] >= 0) {
      sprintf(buf, "%4.2f", prob[i]);
      fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf + 1);
    } else {
      fprintf(fplog, "     ");
    }
  }
  fprintf(fplog, "\n");
}

void ParallelTemperingUtilities::print_count(FILE *fplog, const char *leg,
                                             int n,
                                             const std::vector<int> &count) {
  int i;

  fprintf(fplog, "Repl %2s ", leg);
  for (i = 1; i < n; i++) {
    fprintf(fplog, " %4d", count[i]);
  }
  fprintf(fplog, "\n");
}

void ParallelTemperingUtilities::print_transition_matrix(
    FILE *fplog, int n, const std::vector<std::vector<int>> &nmoves,
    const std::vector<int> &nattempt) {
  int i, j, ntot;
  float Tprint;

  ntot = nattempt[0] + nattempt[1];
  fprintf(fplog, "\n");
  fprintf(fplog, "Repl");
  for (i = 0; i < n; i++) {
    fprintf(fplog, "    "); /* put the title closer to the center */
  }
  fprintf(fplog, "Empirical Transition Matrix\n");

  fprintf(fplog, "Repl");
  for (i = 0; i < n; i++) {
    fprintf(fplog, "%8d", (i + 1));
  }
  fprintf(fplog, "\n");

  for (i = 0; i < n; i++) {
    fprintf(fplog, "Repl");
    for (j = 0; j < n; j++) {
      Tprint = 0.0;
      if (nmoves[i][j] > 0) {
        Tprint = nmoves[i][j] / (2.0 * ntot);
      }
      fprintf(fplog, "%8.4f", Tprint);
    }
    fprintf(fplog, "%3d\n", i);
  }
}

void ParallelTemperingUtilities::print_replica_exchange_statistics(
    FILE *fplog) {
  int i;
  std::vector<bool> nullVec;

  fprintf(fplog, "\nReplica exchange statistics\n");
  if (!bMultiEx) {
    fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
            nattempt[0] + nattempt[1], nattempt[1], nattempt[0]);

    fprintf(fplog, "Repl  average probabilities:\n");
    for (i = 1; i < ms->worldSize; i++) {
      if (nattempt[i % 2] == 0) {
        exchangeProbabilities[i] = 0;
      } else {
        exchangeProbabilities[i] = prob_sum[i] / nattempt[i % 2];
      }
    }
    print_ind(fplog, "", ms->worldSize, ind, nullVec);
    print_prob(fplog, "", ms->worldSize, exchangeProbabilities);

    fprintf(fplog, "Repl  number of exchanges:\n");
    print_ind(fplog, "", ms->worldSize, ind, nullVec);
    print_count(fplog, "", ms->worldSize, nexchange);

    fprintf(fplog, "Repl  average number of exchanges:\n");
    for (i = 1; i < ms->worldSize; i++) {
      if (nattempt[i % 2] == 0) {
        exchangeProbabilities[i] = 0;
      } else {
        exchangeProbabilities[i] =
            (static_cast<double>(nexchange[i])) / nattempt[i % 2];
      }
    }
    print_ind(fplog, "", ms->worldSize, ind, nullVec);
    print_prob(fplog, "", ms->worldSize, exchangeProbabilities);

    fprintf(fplog, "\n");
  }
  /* print the transition matrix */
  print_transition_matrix(fplog, ms->worldSize, nmoves, nattempt);
}

void ParallelTemperingUtilities::print_allswitchind(
    FILE *fplog, int n, const std::vector<int> &pind,
    std::vector<int> &allswaps, std::vector<int> &tmpswap) {
  int i;

  for (i = 0; i < n; i++) {
    tmpswap[i] = allswaps[i];
  }
  for (i = 0; i < n; i++) {
    allswaps[i] = tmpswap[pind[i]];
  }

  fprintf(fplog, "\nAccepted Exchanges:   ");
  for (i = 0; i < n; i++) {
    fprintf(fplog, "%d ", pind[i]);
  }
  fprintf(fplog, "\n");

  /* the "Order After Exchange" is the state label corresponding to the
     configuration that started in state listed in order, i.e.

     3 0 1 2

     means that the:
     configuration starting in simulation 3 is now in simulation 0,
     configuration starting in simulation 0 is now in simulation 1,
     configuration starting in simulation 1 is now in simulation 2,
     configuration starting in simulation 2 is now in simulation 3
   */
  fprintf(fplog, "Order After Exchange: ");
  for (i = 0; i < n; i++) {
    fprintf(fplog, "%d ", allswaps[i]);
  }
  fprintf(fplog, "\n\n");
}

#endif