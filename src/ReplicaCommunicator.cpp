#include "ReplicaCommunicator.h"

#if GOMC_LIB_MPI
ReplicaCommunicator::ReplicaCommunicator(){}

void ReplicaCommunicator::exchangeXYZArrayNonBlocking(XYZArray * myXYZArray, int exchangePartner)
{
    XYZArray outBuffer(*myXYZArray);
    MPI_Request mpi_req;

    int myXYZArrayCount = outBuffer.Count();
    int otherXYZCount;

    MPI_Isend(&myXYZArrayCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&otherXYZCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    XYZArray inBuffer(otherXYZCount);

    MPI_Isend(outBuffer.x, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.x, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.y, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.y, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.z, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.z, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    swap(*myXYZArray, inBuffer);

}

#endif