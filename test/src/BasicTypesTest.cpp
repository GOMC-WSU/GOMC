#include "BasicTypes.h"
#include <gtest/gtest.h>

TEST(BasicTypesTest, CheckDefaultConstructorTest) {
  XYZ point;
  EXPECT_EQ(point.x, 0.0);
  EXPECT_EQ(point.y, 0.0);
  EXPECT_EQ(point.z, 0.0);
}

TEST(BasicTypesTest, CheckConstructorTest) {
  XYZ point(1.0, 99.012, -991.594);
  EXPECT_EQ(point.x, 1.0);
  EXPECT_EQ(point.y, 99.012);
  EXPECT_EQ(point.z, -991.594);
}

TEST(BasicTypesTest, EqualOperatorTest) {
  XYZ point1(1.0, 99.012, -991.594);
  XYZ point2 = point1;
  EXPECT_EQ(point2.x, 1.0);
  EXPECT_EQ(point2.y, 99.012);
  EXPECT_EQ(point2.z, -991.594);
}

TEST(BasicTypesTest, NotEqualOperatorTest) {
  XYZ point1(1.0, 231.2014, 120.18);
  XYZ point2(1.0, 231.2014, 120.18);
  XYZ point3(1.3, 9123.12, 12.18);
  XYZ point4(1.3, 9123.120, 12.18);
  XYZ point5(1.3, 9123.1201, 12.18);

  EXPECT_FALSE(point1 != point2);
  EXPECT_TRUE(point1 != point3);
  EXPECT_FALSE(point3 != point4);
  EXPECT_TRUE(point3 != point5);
}
