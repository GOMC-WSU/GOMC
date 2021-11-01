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
TEST(CheckpointTest, CheckAR_KR) {
    chdir("./test/input/Systems/AR_KR/Base/");
    Simulation base("in.conf");
    base.RunSimulation();

    chdir("../K_1");
    Simulation K_1("in.conf");
    K_1.RunSimulation();  
    chdir("../K_N");
    Simulation K_N("in.conf");

    //K_N.RunSimulation();
    //chdir("../../../../..");
    
    MoleculeLookup & base_ml = base.GetMolLookup();
    for (int i = 0; i < base_ml.molLookupCount; ++i){
        EXPECT_EQ(base_ml.restartMoleculeIndices[i] == base_ml.molLookup[i], true);
        EXPECT_EQ(base_ml.permutedMoleculeIndices[i] == base_ml.molLookup[i], true);
    }
    for (int i = 0; i < base_ml.molLookupCount; ++i){
        EXPECT_EQ(base_ml.restartMoleculeIndices[i] == base_ml.molLookup[i], true);
        EXPECT_EQ(base_ml.permutedMoleculeIndices[i] == base_ml.molLookup[i], true);
    }
        for (int i = 0; i < base_ml.molLookupCount; ++i){
        EXPECT_EQ(base_ml.restartMoleculeIndices[i] == base_ml.molLookup[i], true);
        EXPECT_EQ(base_ml.permutedMoleculeIndices[i] == base_ml.molLookup[i], true);
    }
    
    for (int i = 0; i < base_ml.molLookupCount; ++i){
        std::cout << base_ml.molLookup[i] << " ";
    }        
    std::cout << std::endl;
    for (int i = 0; i < base_ml.molLookupCount; ++i){
        std::cout << base_ml.permutedMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < base_ml.molLookupCount; ++i){
        std::cout << base_ml.restartMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;

    MoleculeLookup & K_1_ml = K_1.GetMolLookup();
    for (int i = 0; i < K_1_ml.molLookupCount; ++i){
        std::cout << K_1_ml.molLookup[i] << " ";
    }        
    std::cout << std::endl;
    for (int i = 0; i < K_1_ml.molLookupCount; ++i){
        std::cout << K_1_ml.permutedMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < K_1_ml.molLookupCount; ++i){
        std::cout << K_1_ml.restartMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;

    MoleculeLookup & K_N_ml = K_N.GetMolLookup();
    for (int i = 0; i < K_N_ml.molLookupCount; ++i){
        std::cout << K_N_ml.molLookup[i] << " ";
    }        
    std::cout << std::endl;
    for (int i = 0; i < K_N_ml.molLookupCount; ++i){
        std::cout << K_N_ml.permutedMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < K_N_ml.molLookupCount; ++i){
        std::cout << K_N_ml.restartMoleculeIndices[i] << " ";
    }
    std::cout << std::endl;
}