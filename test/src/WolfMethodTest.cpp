#include <gtest/gtest.h>
#include "Simulation.h"
#include "FFWolf.h"
#include "NumLib.h"
TEST(WolfMethodTest, CheckLennardJones) {
    EXPECT_EQ(true, true);
}

TEST(WolfMethodTest, CheckElectrostatic) {
    Simulation sim1("test/input/Wolf/in.conf");
    
    // Sets equidistant points with double precision
    sim1.GetCoordinates().Set(0, 0.0, 0.0, 0.0);
    sim1.GetCoordinates().Set(1, 0.0, 1.0, 0.0);
    sim1.GetCoordinates().Set(2, 0.5, sqrt(3)/2.0, 0.0);
    sim1.GetCoordinates().Set(3, 0.5, sqrt(3)/4.0, 3.0/4.0);

    StaticVals * sv = sim1.GetStaticVals();

    double energy = 0.0;

    double r_ij_0_1 = sim1.GetCoordinates().Norm(0,1);
    double r_ij_0_2 = sim1.GetCoordinates().Norm(0,2);
    double r_ij_0_3 = sim1.GetCoordinates().Norm(0,3);
    double r_ij_1_2 = sim1.GetCoordinates().Norm(1,2);
    double r_ij_1_3 = sim1.GetCoordinates().Norm(1,3);
    double r_ij_2_3 = sim1.GetCoordinates().Norm(2,3);

    double wolfAlpha_box_0 = sv->forcefield.wolfAlpha[0];
    double Rc_box_0 = sv->forcefield.rCutCoulomb[0];

    double q_0 = sim1.GetCalcEn().GetCharge(0);
    double q_1 = sim1.GetCalcEn().GetCharge(1);
    double q_2 = sim1.GetCalcEn().GetCharge(2);
    double q_3 = sim1.GetCalcEn().GetCharge(3);

    double summand_0_1 = q_0*q_1*((erfc(wolfAlpha_box_0*r_ij_0_1)/r_ij_0_1) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));
    double summand_0_2 = q_0*q_2*((erfc(wolfAlpha_box_0*r_ij_0_2)/r_ij_0_2) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));
    double summand_0_3 = q_0*q_3*((erfc(wolfAlpha_box_0*r_ij_0_3)/r_ij_0_3) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));
    double summand_1_2 = q_1*q_2*((erfc(wolfAlpha_box_0*r_ij_1_2)/r_ij_1_2) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));
    double summand_1_3 = q_1*q_3*((erfc(wolfAlpha_box_0*r_ij_1_3)/r_ij_1_3) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));
    double summand_2_3 = q_2*q_3*((erfc(wolfAlpha_box_0*r_ij_2_3)/r_ij_2_3) - (erfc(wolfAlpha_box_0*Rc_box_0)/Rc_box_0));

    double sum = summand_0_1 + summand_0_2 + summand_0_3 + summand_1_2 + summand_1_3 + summand_2_3;
    sum *= num::qqFact;

    SystemPotential sysPot = sim1.GetCalcEn().SystemTotal();

    EXPECT_DOUBLE_EQ(sum, sysPot.totalEnergy.totalElect);
}

TEST(WolfMethodTest, CheckLRC) {
    EXPECT_EQ(true, true);
}

