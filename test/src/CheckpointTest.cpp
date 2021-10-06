#include <gtest/gtest.h>
#include "PDBSetup.h"
#include "MolSetup.h"
#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFSetup.h"        //For geometry kinds
#include "FFConst.h"
#include "Reader.h"
#include "InputFileReader.h"
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
TEST(CheckpointTest, CheckAR_KR) {
    chdir("./test/input/Checkpoint/AR_KR/Base/");
    Simulation base("in.conf");
    base.RunSimulation();
    chdir("../K_1");
    Simulation K_1("in.conf");
    K_1.RunSimulation();
    chdir("../K_N");
    Simulation K_N("in.conf");
    K_N.RunSimulation();
    chdir("../../../../..");
}