#include <gtest/gtest.h>
#include "Endian.h"

TEST(EndianTest, TestBitSwap) {
  uint64_t x = 0x0123456789abcdefull;
  uint64_t rev_x = 0xfedcba9876543210ull;
  EXPECT_EQ(rev_x, bswap_64(x));

  uint32_t y = 0x01234567;
  uint32_t rev_y = 0x76543210;
  EXPECT_EQ(rev_y, bswap_32(y));

  uint16_t z = 0x0123;
  uint16_t rev_z = 0x3210;
  EXPECT_EQ(rev_z, bswap_16(z));
}