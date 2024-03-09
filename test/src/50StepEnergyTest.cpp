#include "Simulation.h"
#include <gtest/gtest.h>
#include <unistd.h>
/* There are 4 cases for restarting from checkpoint

1) Base Case:
    The first run.
    Checkpoint Setup  - RestartFromCheckpoint = false
    Checkpoint Output - RestartFromCheckpoint = false

2) The first restart from checkpoint
    Checkpoint Setup  - RestartFromCheckpoint = true
    Checkpoint Output - RestartFromCheckpoint = true

3) The Nth restart from checkpoint from a prev checkpoint
    Checkpoint Setup  - RestartFromCheckpoint = true
    Checkpoint Output - RestartFromCheckpoint = true

*/
TEST(50StepEnergyTest, CheckAR_KR) {
  chdir("./test/input/Systems/AR_KR/Base/");
  Simulation base("in.conf");
  base.RunSimulation();
  /*
  chdir("../K_1");
  Simulation K_1("in.conf");
  K_1.RunSimulation();
  chdir("../K_N");
  Simulation K_N("in.conf");
  K_N.RunSimulation();
  chdir("../../../../..");
  */
  double total = base.GetSystemEnergy();
  // Run the main branch once
  double x2 = 123.00;
  EXPECT_EQ(total, x2);
}