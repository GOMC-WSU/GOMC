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

    chdir("./test/input/Systems/PEN_HEX/Base/");
    Simulation base("in.conf");
    base_runsteps = base.GetRunSteps();
    base.RunSimulation();
    chdir("../Continued");
    Simulation Continued("in.conf");
    chdir("../SingleRun");
    Simulation SingleRun("in100.conf");
    SingleRun.RunSimulation();

    MoleculeLookup & Continued_ml = Continued.GetMolLookup();
    MoleculeLookup & SingleRun_ml = SingleRun.GetMolLookup();
    MoveSettings & Continued_ms = Continued.GetMoveSettings();
    MoveSettings & SingleRun_ms = SingleRun.GetMoveSettings();
    EXPECT_EQ(Continued_ml == SingleRun_ml, true);
    EXPECT_EQ(Continued_ms == SingleRun_ms, true);
    chdir("../../../../..");
    chdir("./test/input/Systems/PEN_HEX");
    system("exec rm -r ./Base/Base_*");
    system("exec rm -r ./Continued/Continued_*");
    system("exec rm -r ./SingleRun/SingleRun_*");
    chdir("../../../..");

}