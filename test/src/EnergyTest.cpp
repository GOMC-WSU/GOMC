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


TEST(EnergyTest, CheckAR_KR) {
    chdir("./test/input/Systems/AR_KR/Base/");
    Simulation base("in.conf");
    base.RunSimulation();
    double total = base.GetSystemEnergy();
    // Run the main branch once
    double x2 = 123.00;
    EXPECT_EQ(total, x2);
}
*/
TEST(EnergyTest, CheckMETHANOL) {
    int result = chdir("./test/input/Systems/METHANOL_OPLSAA/Standard");
    Simulation base("in_GCMC.conf");
    base.RunSimulation();
    double total = base.GetSystemEnergy().Total();
    // Run the main branch once
    double x2 = 6.5846e+06;
    double tol = 0.0001e+06;

    EXPECT_EQ(total, x2);
    EXPECT_TRUE((total >= x2 - tol) && (total <= x2 + tol)); // a is between 1 and 3 inclusive
}