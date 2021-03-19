#include <gtest/gtest.h>
#include "GOMC_Config.h"    //For PT
#include "ReplicaCommunicator.h"
#include "Simulation.h"
#include "ParallelTemperingPreprocessor.h"
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

TEST(ParallelTemperingTest, Pos_And_COMCommunication) {  /// Then you can create tests as usual,
  //using namespace mpi;
  //ompi_communicator_t world;  /// and use MPI inside your tests.

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
    replcomm.exchangeXYZArrayNonBlocking(&newCoords, 1);
    replcomm.exchangeXYZArrayNonBlocking(&newComs, 1);
  } else if (worldRank == 1) {
    replcomm.exchangeXYZArrayNonBlocking(&newCoords, 0);
    replcomm.exchangeXYZArrayNonBlocking(&newComs, 0);
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
    replcomm.exchangeXYZArrayNonBlocking(&newCoords, 1);
    replcomm.exchangeXYZArrayNonBlocking(&newComs, 1);
  } else if (worldRank == 1) {
    replcomm.exchangeXYZArrayNonBlocking(&newCoords, 0);
    replcomm.exchangeXYZArrayNonBlocking(&newComs, 0);
  }

  EXPECT_EQ(oldCoords, newCoords);
  EXPECT_EQ(oldComs, newComs);

}

TEST(ParallelTemperingTest, FullSwapNoEwald) {  /// Then you can create tests as usual,
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

  /* Dummy object, ms->worldRank is corrupted after first swap, hence why we pass it as an argument to ExRep */
  const MultiSim ms(worldSize, worldRank);

  Simulation * sim;

  if(worldRank == 0){
    sim = new Simulation("test/input/ParallelTempering/temp_120.00/repl0.conf", &ms);
  } else if(worldRank == 1){
    sim = new Simulation("test/input/ParallelTempering/temp_180.00/repl1.conf", &ms);
  } else {
    std::cout << worldRank << "something weird happened. " << std::endl;
  }

  Coordinates originalCoords = sim->getCoordinates();
  COM originalCOM = sim->getCOMs();
  //CellList originalCellList = sim->getCellList();
  double originalEnergy = sim->GetSystemEnergy();

  sim->ExchangeReplicas(worldRank);

  Coordinates otherCoords = sim->getCoordinates();
  COM otherCOM = sim->getCOMs();
  //CellList otherCellList = sim->getCellList();
  double otherEnergy = sim->GetSystemEnergy();

  ASSERT_NE(originalCoords, otherCoords);
  ASSERT_NE(originalCOM, otherCOM);
  //ASSERT_NE(originalCellList, otherCellList);
  ASSERT_NE(originalEnergy, otherEnergy);
  
  sim->ExchangeReplicas(worldRank);  

  Coordinates shouldBeOriginalCoords = sim->getCoordinates();
  COM shouldBeOriginalCOM = sim->getCOMs();
  //CellList shouldBeOriginalCellList = sim->getCellList();
  double shouldBeOriginalEnergy = sim->GetSystemEnergy();

  EXPECT_EQ(originalCoords, shouldBeOriginalCoords);
  EXPECT_EQ(originalCOM, shouldBeOriginalCOM);
  //EXPECT_EQ(originalCellList, shouldBeOriginalCellList);
  EXPECT_EQ(originalEnergy, shouldBeOriginalEnergy);

}

#endif

