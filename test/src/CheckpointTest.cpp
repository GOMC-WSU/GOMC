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


TEST(CheckpointTest, CheckPEN_HEX) {

    ulong base_runsteps, Continued_runsteps;
    ulong Continued_true_step;

    int result = chdir("./test/input/Systems/PEN_HEX/Base/");
    Simulation base("in.conf");
    base_runsteps = base.GetRunSteps();
    base.RunSimulation();
    result = chdir("../Continued");
    Simulation Continued("in.conf");
    result = chdir("../SingleRun");
    Simulation SingleRun("in100.conf");
    SingleRun.RunSimulation();

    MoleculeLookup & Continued_ml = Continued.GetMolLookup();
    MoleculeLookup & SingleRun_ml = SingleRun.GetMolLookup();
    MoveSettings & Continued_ms = Continued.GetMoveSettings();
    MoveSettings & SingleRun_ms = SingleRun.GetMoveSettings();

    EXPECT_EQ(Continued_ml.operator==(SingleRun_ml), true);
    EXPECT_EQ(Continued_ms.operator==(SingleRun_ms), true);

    result = chdir("../../../../..");
    result = chdir("./test/input/Systems/PEN_HEX");
    result = system("exec rm -r ./Base/Base_*");
    result = system("exec rm -r ./Continued/Continued_*");
    result = system("exec rm -r ./SingleRun/SingleRun_*");
    result = chdir("../../../..");

}