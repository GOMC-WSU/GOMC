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
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ParallelTemperingUtilities_H
#define ParallelTemperingUtilities_H

#include "GOMC_Config.h"    //For version number
#if GOMC_LIB_MPI
#include <mpi.h>
#endif
#include "XYZArray.h"
#include "ParallelTemperingPreprocessor.h"
#include "System.h"
#include "ConfigSetup.h"

class ParallelTemperingUtilities
{
public:

#if GOMC_LIB_MPI
explicit ParallelTemperingUtilities(MultiSim const*const& multisim, System & sys, StaticVals const& statV, ulong parallelTempFreq, ulong parallelTemperingAttemptsPerExchange);
void evaluateExchangeCriteria(ulong step);
void prepareToDoExchange(const int replica_id, int* maxswap, bool* bThisReplicaExchanged);
void cyclicDecomposition(const vector<int> destinations, vector< vector<int> > & cyclic, vector<bool> & incycle, const int nrepl, int * nswap);
void computeExchangeOrder(vector< vector<int> > & cyclic, vector< vector<int> > & order, const int nrepl, const int maxswap);

void conductExchanges(Coordinates & coordCurrRef, COM & comCurrRef, MultiSim const*const& ms, const int & maxSwap, const bool & bThisReplicaExchanged);
void exchangePositionsNonBlocking(Coordinates & myPos, MultiSim const*const& multisim, int exchangePartner);
void exchangeCOMsNonBlocking(COM & myCOMs, MultiSim const*const& multisim, int exchangePartner);

void exchangePositions(Coordinates & myPos, MultiSim const*const& multisim, int exchangePartner, bool leader);
void exchangeCOMs(COM & myCOMs, MultiSim const*const& multisim, int exchangePartner, bool leader);
void exchangeCellLists(CellList & myCellList, MultiSim const*const& multisim, int exchangePartner, bool leader);
void exchangePotentials(SystemPotential & mySystemPotential, MultiSim const*const& multisim, int exchangePartner, bool leader);
void exchangeVirials(SystemPotential & mySystemPotential, MultiSim const*const& multisim, int exchangePartner, bool leader);
void print_ind(FILE * fplog, const char* leg, int n, vector<int> ind, const vector<bool> bEx);
void print_prob(FILE * fplog, const char* leg, int n, vector<double> prob);
void print_count(FILE * fplog, const char* leg, int n, vector<int> count);
void print_transition_matrix(FILE * fplog, int n, vector< vector<int> > nmoves, const vector<int> nattempt);
void print_replica_exchange_statistics(FILE * fplog);
void print_allswitchind(FILE* fplog, int n, vector<int> pind, vector<int> allswaps, vector<int> tmpswap);
#endif

private:

MultiSim const*const& ms;
bool bMultiEx;
FILE * fplog;
PRNG & prng;
SystemPotential & sysPotRef;
SystemPotential sysPotNew;
ulong parallelTempFreq, parallelTemperingAttemptsPerExchange;
vector<double> global_betas;
vector<int> ind, pind;
vector<bool> exchangeResults;
vector<double> exchangeProbabilities;
vector<int> nexchange;
vector<int> nattempt;

/* For multiex */
vector<int> allswaps;
vector<int> tempswap;
vector< vector<int> > cyclic;
vector< vector<int> > order;
vector<bool> incycle;
/* For multiex */



//! Sum of probabilities
vector<double> prob_sum;
//! Number of moves between replicas i and j
vector< vector<int> > nmoves;

Coordinates newMolsPos;
COM newCOMs;

#if BOX_TOTAL == 1    
    vector<double> global_energies;
#else
    vector<vector <double> > global_energies;
#endif

};

#endif /*ParallelTemperingUtilities_H*/
