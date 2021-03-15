#include <gtest/gtest.h>
#include "GOMC_Config.h"    //For PT
#include "ReplicaCommunicator.h"
#include <iostream>
#include <stdlib.h>
#if GOMC_LIB_MPI
#include <mpi.h>
#include "gtest-mpi-listener.hpp"

GTEST_API_ int main(int argc, char **argv) {
  printf("Running main() from %s\n", __FILE__);
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);

  // Add object that will finalize MPI on exit; Google Test owns this pointer
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);

  // Get the event listener list.
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  // Remove default listener: the default printer and the default XML printer
  ::testing::TestEventListener *l =
        listeners.Release(listeners.default_result_printer());

  // Adds MPI listener; Google Test owns this pointer
  listeners.Append(
      new GTestMPIListener::MPIWrapperPrinter(l,
                                              MPI_COMM_WORLD)
      );
  // Run tests, then clean up and exit. RUN_ALL_TESTS() returns 0 if all tests
  // pass and 1 if some test fails.
  int result = RUN_ALL_TESTS();

  return 0;  // Run tests, then clean up and exit

}

TEST(ParallelTemperingTest, ParallelTemperingTest) {  /// Then you can create tests as usual,
  //using namespace mpi;
  //ompi_communicator_t world;  /// and use MPI inside your tests.
  /* ... test stuff here ... */

  // Get the number of processes
  int worldSize;
  // Initialize the MPI environment
  int worldRank;
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

  std::cout << worldRank << std::endl;

  XYZArray oldCoords(5);
  XYZArray oldComs(5);

  for (int i = 0; i < 5; i++){
    oldCoords.Set(i, (double)(worldRank), (double)(worldRank), (double)(worldRank));
    oldComs.Set(i, (double)(worldRank), (double)(worldRank), (double)(worldRank));
  }

  XYZArray newCoords(5);
  XYZArray newComs(5);

  newCoords = oldCoords;
  newComs = oldComs;

  EXPECT_EQ(oldCoords, newCoords);
  EXPECT_EQ(oldComs, newComs);

  ReplicaCommunicator replcomm;

  if(worldRank == 0){
    for (int i = 0; i < 5; i++)
      std::cout << newCoords[i];
  }

  if(worldRank == 0){
    replcomm.exchangePositionsNonBlocking(&newCoords, 1);
    replcomm.exchangeCOMsNonBlocking(&newComs, 1);
  } else if (worldRank == 1) {
    replcomm.exchangePositionsNonBlocking(&newCoords, 0);
    replcomm.exchangeCOMsNonBlocking(&newComs, 0);
  }

  if(worldRank == 0){
    for (int i = 0; i < 5; i++)
      std::cout << newCoords[i];
  }

  XYZArray interCoords(5);
  XYZArray interComs(5);

  interCoords = newCoords;
  interComs = newComs;

  ASSERT_NE(oldCoords, interCoords);
  ASSERT_NE(oldComs, interComs);

  if(worldRank == 0){
    replcomm.exchangePositionsNonBlocking(&newCoords, 1);
    replcomm.exchangeCOMsNonBlocking(&newComs, 1);
  } else if (worldRank == 1) {
    replcomm.exchangePositionsNonBlocking(&newCoords, 0);
    replcomm.exchangeCOMsNonBlocking(&newComs, 0);
  }

  EXPECT_EQ(oldCoords, newCoords);
  EXPECT_EQ(oldComs, newComs);

}

#endif

