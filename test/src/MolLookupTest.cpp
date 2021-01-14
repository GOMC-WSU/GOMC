#include <gtest/gtest.h>
#include "MoleculeLookup.h"
#include "MoleculeLookup.cpp"


TEST(MolLookupTest, CheckConsensusBeta) {
    /*
        __betas__ molBeta
        1   0   2 -> 1
        2   1   0 -> 1
        2   0   2 -> 2
        0   0   0 -> 0
        2   2   2 -> 2
    */
    uint molBeta;
    MoleculeLookup ml;
    std::vector<double> beta{ 1.0, 0.0, 2.0 };
    molBeta = ml.GetConsensusMolBeta(0, 2, beta, 0, 0, "test");
    EXPECT_EQ(molBeta, 1);

    beta.clear();

    beta.push_back(2.0);
    beta.push_back(1.0);
    beta.push_back(0.0);
    molBeta = ml.GetConsensusMolBeta(0, 2, beta, 0, 0, "test");
    EXPECT_EQ(molBeta, 1);

    beta.clear();

    beta.push_back(2.0);
    beta.push_back(0.0);
    beta.push_back(2.0);
    molBeta = ml.GetConsensusMolBeta(0, 2, beta, 0, 0, "test");
    EXPECT_EQ(molBeta, 2);
    
    beta.clear();

    beta.push_back(0.0);
    beta.push_back(0.0);
    beta.push_back(0.0);
    molBeta = ml.GetConsensusMolBeta(0, 2, beta, 0, 0, "test");
    EXPECT_EQ(molBeta, 0);
        
    beta.clear();

    beta.push_back(2.0);
    beta.push_back(2.0);
    beta.push_back(2.0);
    molBeta = ml.GetConsensusMolBeta(0, 2, beta, 0, 0, "test");
    EXPECT_EQ(molBeta, 2);
}
