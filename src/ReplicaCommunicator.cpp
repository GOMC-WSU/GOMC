#include "ReplicaCommunicator.h"

ReplicaCommunicator::ReplicaCommunicator(){}

void ReplicaCommunicator::exchangePositionsNonBlocking(XYZArray * myPos, int exchangePartner)
{

  XYZArray buffer(*myPos);

  MPI_Request mpi_req;

  MPI_Isend(buffer.x, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->x, myPos->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.y, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->y, myPos->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.z, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->z, myPos->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

}

void ReplicaCommunicator::exchangeCOMsNonBlocking(XYZArray * myCOMs, int exchangePartner)
{

  XYZArray buffer(*myCOMs);

  MPI_Request mpi_req;

  MPI_Isend(buffer.x, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->x, myCOMs->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.y, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->y, myCOMs->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer.z, buffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->z, myCOMs->Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

}