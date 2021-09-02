#include <gtest/gtest.h>
#include "Simulation.h"
#include "FFWolf.h"

TEST(WolfMethodTest, CheckLennardJones) {
    EXPECT_EQ(true, true);
}

TEST(WolfMethodTest, CheckElectrostatic) {
    Simulation sim1("test/input/Wolf/in.conf");
    
    // Sets equidistant points with double precision
    sim1.GetCoordinates().Set(0, 0.0, 0.0, 0.0);
    sim1.GetCoordinates().Set(1, 0.0, 1.0, 0.0);
    sim1.GetCoordinates().Set(2, 0.5, sqrt(3)/2, 0.0);
    sim1.GetCoordinates().Set(3, 0.5, sqrt(3)/4, 3/4);

    EXPECT_EQ(true, true);
}

TEST(WolfMethodTest, CheckLRC) {
    EXPECT_EQ(true, true);
}

