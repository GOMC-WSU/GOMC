#include "BitLib.h"
#include <climits>
#include <gtest/gtest.h>

TEST(BitLabTest, CheckTest) {
  uint v = 8;
  uint pos = 3;
  EXPECT_EQ(bits::Check(v, pos), v);
  EXPECT_EQ(bits::Check(v, pos - 1), 0);
  EXPECT_EQ(bits::Check(v, pos + 1), 0);
  EXPECT_EQ(bits::Check(UINT_MAX, 30), 0x40000000);
  EXPECT_EQ(bits::Check(UINT_MAX, 28), 0x10000000);
}

TEST(BitLabTest, CountSetTest) {
  uint v = 336070675;
  EXPECT_EQ(bits::CountSet(v), 7);
  EXPECT_EQ(bits::CountSet(UINT_MAX), 32);
}

TEST(BitLabTest, GetMasksTest) {
  uint v = 4;
  std::vector<std::vector<uint>> mask{
      {0}, {1, 2, 4, 8}, {3, 5, 6, 9, 10, 12}, {7, 11, 13, 14}};
  EXPECT_EQ(bits::GetMasks(v), mask);

  v = 0;
  mask.clear();
  EXPECT_EQ(bits::GetMasks(0), mask);

  std::vector<std::vector<uint>> mask1{{0}};
  EXPECT_EQ(bits::GetMasks(1), mask1);

  std::vector<std::vector<uint>> mask2{{0}, {1, 2}};
  EXPECT_EQ(bits::GetMasks(2), mask2);
}
