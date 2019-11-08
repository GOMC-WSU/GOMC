/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PARALLELTEMPERINGMPIMETHODS_H
#define PARALLELTEMPERINGMPIMETHODS_H
#include <mpi.h>
#include "GOMC_Config.h"

class ParallelTemperingMPIMethods{
//    ParallelTemperingMPIMethods();
public:
static void gomc_sumd_comm(int nr, double r[], MPI_Comm mpi_comm);
};
#endif