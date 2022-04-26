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

    // This is needed to get passed Remark 
    uint frameNum = 1;
    pdbBase.Init(rsBase, pdbnamesBase, frameNum);    
    pdbRefBase.Init(rsBaseRef, pdbnamesBaseRef, frameNum);

    frameNum = 11;
    pdbBase.Init(rsBase, pdbnamesBase, frameNum);    
    pdbBaseRef.Init(rsBaseRef, pdbnamesBaseRef, frameNum);

    config_setup::RestartSettings rsBaseRestart;
    config_setup::RestartSettings rsBaseRefRestart;

    rsBaseRestart.restartFromCheckpoint = true;
    rsBaseRestartRef.restartFromCheckpoint = true;

    pdbBaseRestart.Init(rsBaseRestart, pdbnamesBaseRestart);    
    pdbBaseRefRestart.Init(rsBaseRefRestart, pdbnamesBaseRefRestart);
    

    // Checks if the coordinates of the traj match the restart file
    for (uint i = 0; i < pdbBase.atoms.count; ++i){
        ptrdiff_t pos = find(pdbBaseRestart.atoms.chainLetter.begin(), 
            pdbBaseRestart.atoms.chainLetter.end(), pdbBase.atoms.chainLetter[i])
            - pdbBaseRestart.atoms.chainLetter.begin();
        if(pdbBase.atoms.x[i] != pdbBaseRestart.atoms.x[pos])
            std::cout << pdbBase.atoms.x[i] << " " << i << " " << pdbBaseRestart.atoms.x[pos] << " " << pos <<  std::endl;
        EXPECT_EQ(pdbBase.atoms.x[i] == pdbBaseRestart.atoms.x[pos], true);
    }

    // Checks if the coordinates of the traj match the restart file
    for (uint i = 0; i < pdb1.atoms.count; ++i){
        ptrdiff_t pos = find(pdb1Restart.atoms.chainLetter.begin(), 
            pdb1Restart.atoms.chainLetter.end(), pdb1.atoms.chainLetter[i])
            - pdb1Restart.atoms.chainLetter.begin();
        if(pdb1.atoms.x[i] != pdb1Restart.atoms.x[pos])
            std::cout << pdbBase.atoms.x[i] << " " << i << " " << pdb1Restart.atoms.x[pos] << " " << pos <<  std::endl;
        EXPECT_EQ(pdb1.atoms.x[i] == pdb1Restart.atoms.x[pos], true);
    }

    // Checks if the coordinates of the traj match the restart file
    for (uint i = 0; i < pdbN.atoms.count; ++i){
        ptrdiff_t pos = find(pdbNRestart.atoms.chainLetter.begin(), 
            pdbNRestart.atoms.chainLetter.end(), pdbN.atoms.chainLetter[i])
            - pdbNRestart.atoms.chainLetter.begin();
        if(pdbN.atoms.x[i] != pdbNRestart.atoms.x[pos])
            std::cout << pdbN.atoms.x[i] << " " << i << " " << pdbNRestart.atoms.x[pos] << " " << pos <<  std::endl;
        EXPECT_EQ(pdbN.atoms.x[i] == pdbNRestart.atoms.x[pos], true);
    }

    PDBSetup pdbBase_Base_To_K_1, pdb1_Base_To_K_1;

    // This is needed to get passed Remark 
    // Not sure why.
    frameNum = 1;
    pdbBase_Base_To_K_1.Init(rsBase, pdbnamesBase, frameNum);    
    pdb1_Base_To_K_1.Init(rs1, pdbnames1, frameNum);

    int lastFrame = 11;
    pdbBase_Base_To_K_1.Init(rsBase, pdbnamesBase, lastFrame);    

    // Checks if the last frame the base traj match the first frame of K_1 traj
    for (uint i = 0; i < pdbBase_Base_To_K_1.atoms.count; ++i){
        // Find mol i's chain index in restart output files 
        ptrdiff_t pos1 = find(pdb1_Base_To_K_1.atoms.chainLetter.begin(), 
            pdb1_Base_To_K_1.atoms.chainLetter.end(), pdbBase_Base_To_K_1.atoms.chainLetter[i])
            - pdb1_Base_To_K_1.atoms.chainLetter.begin();
         ptrdiff_t pos2 = find(pdbBase_Base_To_K_1.atoms.chainLetter.begin(), 
            pdbBase_Base_To_K_1.atoms.chainLetter.end(), pdb1_Base_To_K_1.atoms.chainLetter[i])
            - pdbBase_Base_To_K_1.atoms.chainLetter.begin();
        if(pdbBase_Base_To_K_1.atoms.x[i] != pdb1_Base_To_K_1.atoms.x[pos1])
            std::cout << pdbBase_Base_To_K_1.atoms.x[i] << " " << i << " " << pdbBaseRestart.atoms.x[pos1] << " " << pos1 <<  std::endl;
        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.x[i] == pdb1_Base_To_K_1.atoms.x[pos1], true);
        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.y[i] == pdb1_Base_To_K_1.atoms.y[pos1], true);
        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.z[i] == pdb1_Base_To_K_1.atoms.z[pos1], true);

        EXPECT_EQ(pdb1_Base_To_K_1.atoms.x[i] == pdbBase_Base_To_K_1.atoms.x[pos2], true);
        EXPECT_EQ(pdb1_Base_To_K_1.atoms.y[i] == pdbBase_Base_To_K_1.atoms.y[pos2], true);
        EXPECT_EQ(pdb1_Base_To_K_1.atoms.z[i] == pdbBase_Base_To_K_1.atoms.z[pos2], true);

        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.chainLetter[i] == pdb1_Base_To_K_1.atoms.chainLetter[pos1], true);
        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.box[i] == pdb1_Base_To_K_1.atoms.box[pos1], true);
        EXPECT_EQ(pdbBase_Base_To_K_1.atoms.resNames[i] == pdb1_Base_To_K_1.atoms.resNames[pos1], true);
        EXPECT_EQ(pdb1_Base_To_K_1.atoms.beta[i] == pdbBase_Base_To_K_1.atoms.beta[pos2], true);

        EXPECT_EQ(pos1 == pos2, true);
    }
    PDBSetup pdb1_K_1_To_K_N, pdnN_K_1_To_K_N;

    // This is needed to get passed Remark 
    // Not sure why.
    frameNum = 1;
    pdb1_K_1_To_K_N.Init(rs1, pdbnames1, frameNum);    
    pdnN_K_1_To_K_N.Init(rsN, pdbnamesN, frameNum);

    lastFrame = 11;
    pdb1_K_1_To_K_N.Init(rs1, pdbnames1, lastFrame);    

    // Checks if the last frame the base traj match the first frame of K_1 traj
    for (uint i = 0; i < pdb1_K_1_To_K_N.atoms.count; ++i){
        // Find mol i's chain index in restart output files 
        ptrdiff_t pos1 = find(pdnN_K_1_To_K_N.atoms.chainLetter.begin(), 
            pdnN_K_1_To_K_N.atoms.chainLetter.end(), pdb1_K_1_To_K_N.atoms.chainLetter[i])
            - pdnN_K_1_To_K_N.atoms.chainLetter.begin();
         ptrdiff_t pos2 = find(pdb1_K_1_To_K_N.atoms.chainLetter.begin(), 
            pdb1_K_1_To_K_N.atoms.chainLetter.end(), pdnN_K_1_To_K_N.atoms.chainLetter[i])
            - pdb1_K_1_To_K_N.atoms.chainLetter.begin();
        if(pdb1_K_1_To_K_N.atoms.x[i] != pdnN_K_1_To_K_N.atoms.x[pos1])
            std::cout << pdb1_K_1_To_K_N.atoms.x[i] << " " << i << " " << pdbBaseRestart.atoms.x[pos1] << " " << pos1 <<  std::endl;
        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.x[i] == pdnN_K_1_To_K_N.atoms.x[pos1], true);
        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.y[i] == pdnN_K_1_To_K_N.atoms.y[pos1], true);
        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.z[i] == pdnN_K_1_To_K_N.atoms.z[pos1], true);

        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.x[i] == pdb1_K_1_To_K_N.atoms.x[pos2], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.y[i] == pdb1_K_1_To_K_N.atoms.y[pos2], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.z[i] == pdb1_K_1_To_K_N.atoms.z[pos2], true);

        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.chainLetter[i] == pdnN_K_1_To_K_N.atoms.chainLetter[pos1], true);
        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.box[i] == pdnN_K_1_To_K_N.atoms.box[pos1], true);
        EXPECT_EQ(pdb1_K_1_To_K_N.atoms.resNames[i] == pdnN_K_1_To_K_N.atoms.resNames[pos1], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.beta[i] == pdb1_K_1_To_K_N.atoms.beta[pos2], true);

        EXPECT_EQ(pos1 == pos2, true);
    }
    
    PDBSetup pdb_SingleRun;
    frameNum = 1;
    pdb_SingleRun.Init(rsSingleRun, pdbnamesSingleRun, frameNum);  
    lastFrame = 11;
    pdnN_K_1_To_K_N.Init(rsN, pdbnamesN, lastFrame);
    lastFrame = 31;
    pdb_SingleRun.Init(rsSingleRun, pdbnamesSingleRun, lastFrame);  

    // Checks if the last frame the SingleRun traj match the last frame of K_N traj
    for (uint i = 0; i < pdb1_K_1_To_K_N.atoms.count; ++i){
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.x[i] == pdb_SingleRun.atoms.x[i], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.y[i] == pdb_SingleRun.atoms.y[i], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.z[i] == pdb_SingleRun.atoms.z[i], true);

        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.chainLetter[i] == pdb_SingleRun.atoms.chainLetter[i], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.box[i] == pdb_SingleRun.atoms.box[i], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.resNames[i] == pdb_SingleRun.atoms.resNames[i], true);
        EXPECT_EQ(pdnN_K_1_To_K_N.atoms.beta[i] == pdb_SingleRun.atoms.beta[i], true);
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

