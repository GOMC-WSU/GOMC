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
    //EXPECT_EQ(Continued_ml == SingleRun_ml, true);
    //EXPECT_EQ(Continued_ms == SingleRun_ms, true);
    EXPECT_EQ(Continued_ms.scale == SingleRun_ms.scale, true);
    EXPECT_EQ(Continued_ms.acceptPercent == SingleRun_ms.acceptPercent, true);
    //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
    //molLookup
    //index [BOX_TOTAL * kind + box + 1] is the element after the end
    //of that kind/box
    EXPECT_EQ(Continued_ms.accepted == SingleRun_ms.accepted, true);
    EXPECT_EQ(Continued_ms.tries == SingleRun_ms.tries, true);
    EXPECT_EQ(Continued_ms.tempAccepted == SingleRun_ms.tempAccepted, true);
    EXPECT_EQ(Continued_ms.tempTries == SingleRun_ms.tempTries, true);
    EXPECT_EQ(Continued_ms.mp_accepted == SingleRun_ms.mp_accepted, true);
    EXPECT_EQ(Continued_ms.mp_tries == SingleRun_ms.mp_tries, true);//Kinds that can move intra and inter box
    EXPECT_EQ(Continued_ms.mp_interval_accepted == SingleRun_ms.mp_interval_accepted, true); //Kinds that can move intra box only
    EXPECT_EQ(Continued_ms.mp_interval_tries == SingleRun_ms.mp_interval_tries, true);// stores the molecule index for global atom index
    EXPECT_EQ(Continued_ms.mp_r_max == SingleRun_ms.mp_r_max, true); // stores the local atom index for global atom index
    EXPECT_EQ(Continued_ms.mp_t_max == SingleRun_ms.mp_t_max, true); // stores the molecule kind for global atom index

    //chdir("../../../../..");

}