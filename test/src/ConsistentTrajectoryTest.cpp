#include <gtest/gtest.h>
#include "PDBSetup.h"
#include "MolSetup.h"
#include "BondAdjacencyList.h"
#include "ConfigSetup.h"
#include "FFSetup.h"        //For geometry kinds
#include "FFConst.h"
#include "Reader.h"
#include "InputFileReader.h"
#include "AlphaNum.h"
#include "CheckpointSetup.h"

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

TEST(ConsistentTrajectoryTest, CheckChainLetters) {

    config_setup::RestartSettings rs;
    config_setup::InFiles if2;

    std::string pdbnamesSTART[2], pdbnamesBase[2], pdbnames1[2], pdbnamesN[2], pdbnamesNPlus1[2];

    pdbnamesSTART[0] = "./test/input/ConsistentTrajectory/START_BOX_0.pdb";
    pdbnamesSTART[1] = "./test/input/ConsistentTrajectory/START_BOX_1.pdb";

    pdbnamesBase[0] = "./test/input/ConsistentTrajectory/base_BOX_0.pdb";
    pdbnamesBase[1] = "./test/input/ConsistentTrajectory/base_BOX_1.pdb";

    pdbnames1[0] = "./test/input/ConsistentTrajectory/1_BOX_0.pdb";
    pdbnames1[1] = "./test/input/ConsistentTrajectory/1_BOX_1.pdb";

    pdbnamesN[0] = "./test/input/ConsistentTrajectory/n_BOX_0.pdb";
    pdbnamesN[1] = "./test/input/ConsistentTrajectory/n_BOX_1.pdb";

    PDBSetup pdbStart, pdbBase, pdb1, pdbN;

    pdbStart.Init(rs, if2, pdbnamesSTART);
    pdbBase.Init(rs, if2, pdbnamesBase);    
    pdb1.Init(rs, if2, pdbnames1);
    pdbN.Init(rs, if2, pdbnamesN);
    
    for (uint i = 0; i < pdbStart.atoms.count; ++i){
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdbBase.atoms.chainLetter[i], true);
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdb1.atoms.chainLetter[i], true);
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdbN.atoms.chainLetter[i], true);
    }
}