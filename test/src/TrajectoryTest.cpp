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

TEST(FiftyStepTrajectoryTest, Check_AR_KR) {

    ulong base_runsteps;
    std::string inputFileString("in.conf");
    int result = chdir("./test/input/Systems/AR_KR/Base/");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    } else {
        Simulation base(inputFileString.c_str());
        base_runsteps = base.GetRunSteps();
        base.RunSimulation();
    }
    result = chdir("../../../../..");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    config_setup::RestartSettings rsStart;
    config_setup::RestartSettings rsBase, rsBaseRef;

    std::string pdbnamesBase[2], pdbnamesBaseRef[2];
    std::string pdbnamesBaseRestart[2], pdbnamesBaseRefRestart[2];

    pdbnamesBase[0] = "./test/input/Systems/AR_KR/Base/base_BOX_0.pdb";
    pdbnamesBase[1] = "./test/input/Systems/AR_KR/Base/base_BOX_1.pdb";

    pdbnamesBaseRef[0] = "./test/input/Systems/AR_KR/BaseRef/base_BOX_0.pdb";
    pdbnamesBaseRef[1] = "./test/input/Systems/AR_KR/BaseRef/base_BOX_0.pdb";

    pdbnamesBaseRestart[0] = "./test/input/Systems/AR_KR/Base/base_BOX_0_restart.pdb";
    pdbnamesBaseRestart[1] = "./test/input/Systems/AR_KR/Base/base_BOX_1_restart.pdb";

    pdbnamesBaseRefRestart[0] = "./test/input/Systems/AR_KR/BaseRef/base_BOX_0_restart.pdb";
    pdbnamesBaseRefRestart[1] = "./test/input/Systems/AR_KR/BaseRef/base_BOX_1_restart.pdb";

    PDBSetup pdbBase, pdbBaseRef;
    PDBSetup pdbBaseRestart, pdbBaseRefRestart;

    rsStart.recalcTrajectory = false;
    rsBase.recalcTrajectory = true;
    rsBaseRef.recalcTrajectory = true;
    
    int frameNum = 1;
    pdbBase.Init(rsBase, pdbnamesBase, frameNum);  
    pdbBaseRef.Init(rsBaseRef, pdbnamesBaseRef, frameNum);
  
    for(int frame = 0; frame < 10; ++frame){
        pdbBaseRef.Init(rsBaseRef, pdbnamesBaseRef, frameNum+frame);
        pdbBase.Init(rsBase, pdbnamesBase, frameNum+frame);    

            // Checks if the last frame the SingleRun traj match the last frame of K_N traj
        for (uint i = 0; i < pdbBaseRef.atoms.count; ++i){
            EXPECT_EQ(pdbBaseRef.atoms.x[i] == pdbBase.atoms.x[i], true);
            EXPECT_EQ(pdbBaseRef.atoms.y[i] == pdbBase.atoms.y[i], true);
            EXPECT_EQ(pdbBaseRef.atoms.z[i] == pdbBase.atoms.z[i], true);

            EXPECT_EQ(pdbBaseRef.atoms.chainLetter[i] == pdbBase.atoms.chainLetter[i], true);
            EXPECT_EQ(pdbBaseRef.atoms.box[i] == pdbBase.atoms.box[i], true);
            EXPECT_EQ(pdbBaseRef.atoms.resNames[i] == pdbBase.atoms.resNames[i], true);
            EXPECT_EQ(pdbBaseRef.atoms.beta[i] == pdbBase.atoms.beta[i], true);
        }
    }
    result = chdir("./test/input/Systems/AR_KR");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = system("exec rm -r ./Base/base_*");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = chdir("../../../..");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
}

