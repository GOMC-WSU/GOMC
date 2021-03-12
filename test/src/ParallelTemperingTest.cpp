#include <gtest/gtest.h>
#include "GOMC_Config.h"    //For PT
#include "Simulation.h"

#if GOMC_LIB_MPI
#include <mpi.h>


TEST(ParallelTemperingTest, ParallelTemperingTest) {  /// Then you can create tests as usual,
  //using namespace mpi;
  //ompi_communicator_t world;  /// and use MPI inside your tests.
  /* ... test stuff here ... */
    EXPECT_EQ(1, 1);
    Simulation sim0("sim0.conf");
    Simulation sim1("sim1.conf");

}

#endif

