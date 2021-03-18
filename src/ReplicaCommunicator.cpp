#include "ReplicaCommunicator.h"

#if GOMC_LIB_MPI
ReplicaCommunicator::ReplicaCommunicator(){}

void ReplicaCommunicator::exchangePositionsNonBlocking(Coordinates * myPos, int exchangePartner)
{

  int count = myPos->Count();


  double * buffer_x = new double[count];
  double * buffer_y = new double[count];
  double * buffer_z = new double[count];
/*
#ifdef _OPENMP
  #pragma omp parallel default(none) 
#endif
  {
    std::memcpy(buffer_x, myPos->x, sizeof(double) * count);
    std::memcpy(buffer_y, myPos->y, sizeof(double) * count);
    std::memcpy(buffer_z, myPos->z, sizeof(double) * count);
  }
*/
    std::memcpy(buffer_x, myPos->x, sizeof(double) * count);
    std::memcpy(buffer_y, myPos->y, sizeof(double) * count);
    std::memcpy(buffer_z, myPos->z, sizeof(double) * count);
  
  MPI_Request mpi_req;

  MPI_Isend(buffer_x, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->x, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer_y, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->y, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer_z, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myPos->z, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  delete buffer_x;
  delete buffer_y;
  delete buffer_z;


}

void ReplicaCommunicator::exchangeCOMsNonBlocking(COM * myCOMs, int exchangePartner)
{
  int count = myCOMs->Count();

  //XYZArray buffer(*myCOMs);
  double * buffer_x = new double[count];
  double * buffer_y = new double[count];
  double * buffer_z = new double[count];
/*
  #ifdef _OPENMP
  #pragma omp parallel default(none) shared(count)
#endif
  {
    std::memcpy(buffer_x, myCOMs->x, sizeof(double) * count);
    std::memcpy(buffer_y, myCOMs->y, sizeof(double) * count);
    std::memcpy(buffer_z, myCOMs->z, sizeof(double) * count);
  }
*/
  
    std::memcpy(buffer_x, myCOMs->x, sizeof(double) * count);
    std::memcpy(buffer_y, myCOMs->y, sizeof(double) * count);
    std::memcpy(buffer_z, myCOMs->z, sizeof(double) * count);
  

 

  MPI_Request mpi_req;

  MPI_Isend(buffer_x, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->x, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer_y, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->y, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  MPI_Isend(buffer_z, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
  MPI_Recv(myCOMs->z, count * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

  delete buffer_x;
  delete buffer_y;
  delete buffer_z;

}

#endif