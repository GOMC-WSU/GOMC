#include "MoleculeLookup.h"
#include <gtest/gtest.h>

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
  // MoleculeLookup ml;
  std::vector<double> beta(3);

  beta[0] = 1.0;
  beta[1] = 0.0;
  beta[2] = 2.0;

  molBeta = MoleculeLookup::GetConsensusMolBeta(0, 2, beta, 0, 0, "test0");
  EXPECT_EQ(molBeta, 1);

  beta[0] = 2.0;
  beta[1] = 1.0;
  beta[2] = 0.0;

  molBeta = MoleculeLookup::GetConsensusMolBeta(0, 2, beta, 0, 0, "test1");
  EXPECT_EQ(molBeta, 1);

  beta[0] = 2.0;
  beta[1] = 0.0;
  beta[2] = 2.0;

  molBeta = MoleculeLookup::GetConsensusMolBeta(0, 2, beta, 0, 0, "test2");
  EXPECT_EQ(molBeta, 2);

  beta[0] = 0.0;
  beta[1] = 0.0;
  beta[2] = 0.0;

  molBeta = MoleculeLookup::GetConsensusMolBeta(0, 2, beta, 0, 0, "test3");
  EXPECT_EQ(molBeta, 0);

  beta[0] = 2.0;
  beta[1] = 2.0;
  beta[2] = 2.0;

  molBeta = MoleculeLookup::GetConsensusMolBeta(0, 2, beta, 0, 0, "test4");
  EXPECT_EQ(molBeta, 2);
}
