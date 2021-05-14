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

    config_setup::RestartSettings rsStart;
    config_setup::RestartSettings rsBase;
    config_setup::RestartSettings rs1;
    config_setup::RestartSettings rsN;

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

    rsStart.recalcTrajectory = false;

    pdbStart.Init(rsStart, if2, pdbnamesSTART);
    pdbBase.Init(rsBase, if2, pdbnamesBase);    
    pdb1.Init(rs1, if2, pdbnames1);
    pdbN.Init(rsN, if2, pdbnamesN);
    
    for (uint i = 0; i < pdbStart.atoms.count; ++i){
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdbBase.atoms.chainLetter[i], true);
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdb1.atoms.chainLetter[i], true);
        EXPECT_EQ(pdbStart.atoms.chainLetter[i] == pdbN.atoms.chainLetter[i], true);
    }
}

TEST(ConsistentTrajectoryTest, CheckPDBTrajCoordinates) {

    config_setup::RestartSettings rsStart;
    config_setup::RestartSettings rsBase;
    config_setup::RestartSettings rs1;
    config_setup::RestartSettings rsN;

    config_setup::InFiles if2;

    std::string pdbnamesSTART[2], pdbnamesBase[2], pdbnames1[2], pdbnamesN[2];
    std::string pdbnamesBaseRestart[2], pdbnames1Restart[2], pdbnamesNRestart[2];

    pdbnamesSTART[0] = "./test/input/ConsistentTrajectory/START_BOX_0.pdb";
    pdbnamesSTART[1] = "./test/input/ConsistentTrajectory/START_BOX_1.pdb";

    pdbnamesBase[0] = "./test/input/ConsistentTrajectory/base_BOX_0.pdb";
    pdbnamesBase[1] = "./test/input/ConsistentTrajectory/base_BOX_1.pdb";

    pdbnames1[0] = "./test/input/ConsistentTrajectory/1_BOX_0.pdb";
    pdbnames1[1] = "./test/input/ConsistentTrajectory/1_BOX_1.pdb";

    pdbnamesN[0] = "./test/input/ConsistentTrajectory/n_BOX_0.pdb";
    pdbnamesN[1] = "./test/input/ConsistentTrajectory/n_BOX_1.pdb";

    pdbnamesBaseRestart[0] = "./test/input/ConsistentTrajectory/base_BOX_0_restart.pdb";
    pdbnamesBaseRestart[1] = "./test/input/ConsistentTrajectory/base_BOX_1_restart.pdb";

    pdbnames1Restart[0] = "./test/input/ConsistentTrajectory/1_BOX_0_restart.pdb";
    pdbnames1Restart[1] = "./test/input/ConsistentTrajectory/1_BOX_1_restart.pdb";

    pdbnamesNRestart[0] = "./test/input/ConsistentTrajectory/n_BOX_0_restart.pdb";
    pdbnamesNRestart[1] = "./test/input/ConsistentTrajectory/n_BOX_1_restart.pdb";

    PDBSetup pdbStart, pdbBase, pdb1, pdbN;
    PDBSetup pdbBaseRestart, pdb1Restart, pdbNRestart;

    rsStart.recalcTrajectory = false;
    rsBase.recalcTrajectory = false;
    rs1.recalcTrajectory = false;
    rsN.recalcTrajectory = false;

    uint frameNum = 1;

    pdbStart.Init(rsStart, if2, pdbnamesSTART);
    pdbBase.Init(rsBase, if2, pdbnamesBase, frameNum);    
    pdb1.Init(rs1, if2, pdbnames1, frameNum);
    pdbN.Init(rsN, if2, pdbnamesN, frameNum);

/*
    config_setup::RestartSettings rsBaseRestart;
    config_setup::RestartSettings rs1Restart;
    config_setup::RestartSettings rsNRestart;

    rsBaseRestart.restartFromCheckpoint;
    rs1Restart.restartFromCheckpoint;
    rsNRestart.restartFromCheckpoint;
*/
    pdbBaseRestart.Init(rsBase, if2, pdbnamesBaseRestart);    
    pdb1Restart.Init(rs1, if2, pdbnames1Restart);
    pdbNRestart.Init(rsN, if2, pdbnamesNRestart);
    
    for (uint i = 0; i < pdbStart.atoms.count; ++i){
        /* Find mol i's chain index in restart output files */
        ptrdiff_t pos = find(pdbBaseRestart.atoms.chainLetter.begin(), 
            pdbBaseRestart.atoms.chainLetter.end(), pdbStart.atoms.chainLetter[i])
            - pdbBaseRestart.atoms.chainLetter.begin();
        if(pdbBase.atoms.x[i] != pdbBaseRestart.atoms.x[pos])
            std::cout << pdbBase.atoms.x[i] << " " << pdbBaseRestart.atoms.x[pos] << std::endl;
        EXPECT_EQ(pdbBase.atoms.x[i] == pdbBaseRestart.atoms.x[pos], true);
    }
}