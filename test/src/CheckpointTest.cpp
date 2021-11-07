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
#include "MoveSettings.h"
#ifdef GOMC_CUDA
#include "cuda.h"
#include <cuda_runtime_api.h>
#endif

TEST(CheckpointTest, CheckPEN_HEX) {
#ifdef GOMC_CUDA
    system("nvidia-smi --gpu-reset");
#endif
    ulong base_runsteps, Continued_runsteps;
    ulong Continued_true_step;

    int result = chdir("./test/input/Systems/PEN_HEX/Base/");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    Simulation base("in.conf");
    base.RunSimulation();
    result = chdir("../Continued");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    Simulation Continued("in.conf");
    result = chdir("../SingleRun");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    Simulation SingleRun("in100.conf");
    SingleRun.RunSimulation();

    MoleculeLookup & Continued_ml = Continued.GetMolLookup();
    MoleculeLookup & SingleRun_ml = SingleRun.GetMolLookup();
    MoveSettings & Continued_ms = Continued.GetMoveSettings();
    MoveSettings & SingleRun_ms = SingleRun.GetMoveSettings();

    EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
    EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);

    result = chdir("../../../../..");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    } 
    result = chdir("./test/input/Systems/PEN_HEX");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = system("exec rm -r ./Base/Base_*");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = system("exec rm -r ./Continued/Continued_*");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = system("exec rm -r ./SingleRun/SingleRun_*");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
    result = chdir("../../../..");
    if (result){
        std::cout << "System call failed!" << std::endl;
        exit(1);
    }
#ifdef GOMC_CUDA
    cudaDeviceReset();
#endif
}