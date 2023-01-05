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
TEST(WolfCalibrationTest, CheckSPCE_ETOH) {
  chdir("./test/input/WolfCalibration");
  Simulation base("wolf_calibration.conf");
  base.RunSimulation(10000);
  double total = base.GetSystemEnergy().totalEnergy.total;
  // Run the main branch once
  double x2 = 123.00;
  EXPECT_EQ(total, x2);
}