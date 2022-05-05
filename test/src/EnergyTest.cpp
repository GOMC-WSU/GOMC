#include <gtest/gtest.h>
#include "Simulation.h"
#include<unistd.h> 
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
    double total = base.GetSystemEnergy();
    // Run the main branch once
    double x2 = 123.00;
    EXPECT_EQ(total, x2);
}

TEST(50StepEnergyTest, CheckMETHANOL) {
    chdir("./test/input/Systems/METHANOL_OPLSAA/Standard");
    Simulation base("in_GCMC.conf");
    base.RunSimulation();
    double total = base.GetSystemEnergy();
    // Run the main branch once
    double x2 = 123.00;
    EXPECT_EQ(total, x2);
}