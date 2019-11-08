/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "ParallelTemperingMPIMethods.h"
#if GOMC_LIB_MPI

void ParallelTemperingMPIMethods::gomc_sumd_comm(int nr, double r[], MPI_Comm mpi_comm)
{
#if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication
       structure), so this isn't as bad as it looks. */
    double *buf;
    int     i;

    buf = new double[nr];
    MPI_Allreduce(r, buf, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
    for (i = 0; i < nr; i++)
    {
        r[i] = buf[i];
    }
    delete[] buf;
#endif
#endif
}