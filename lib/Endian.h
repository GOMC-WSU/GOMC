#pragma once

/*
 * This file was created to handle endianness problem when writing to binary
 * files. The bug could occur when an integer is being written to a binary file
 * and read from a different system. If systems have different endianness, the
 * order of bytes could get reversed.
 *
 * NUM: 0x04030201
 * LITTLE_ENDIAN => 0x01020304
 * BIG_ENDIAN    => 0x04030201
 *
 * e.g. If the writing system (little endian) have 0x04030201 in its variable
 * and wants to write it to file, it would be in the order of 0x01, 0x02, 0x03
 * 0x04. And when the reading system (big endian) reads with the same order and
 * store it in an integer, our original 0x04030201 becomes 0x01020304.
 *
 * To prevent this from happening, we are always going to assume that the file
 * format is going to be little endian, and GOMC will perform conversion to
 * little or big endian based on the type of system.
 *
 * This header file will include tools to detect the endianness of the system
 * and functions to convert the little endian to big endian and big endian to
 * little endian.
 *
 * For developers outside of this header file, two functions should only be
 * used:
 *
 * // converts host integer to file integer
 * uint64_t htof64(uint64_t host_integer)
 * // converts host integer to file integer
 * uint64_t ftoh64(uint64_t file_integer)
 *
 * Before writing to file, make sure you use htof64() and after reading use
 * ftoh64()!
 *
 */

#include <stdlib.h>
#include <stdint.h>

#define bswap_64(x)                                                            \
  ((((x)&0xff00000000000000ull) >> 56) | (((x)&0x00ff000000000000ull) >> 40) | \
   (((x)&0x0000ff0000000000ull) >> 24) | (((x)&0x000000ff00000000ull) >> 8) |  \
   (((x)&0x00000000ff000000ull) << 8) | (((x)&0x0000000000ff0000ull) << 24) |  \
   (((x)&0x000000000000ff00ull) << 40) | (((x)&0x00000000000000ffull) << 56))

#define bswap_32(x)                                                            \
  ((((x)&0xff000000) >> 24) | (((x)&0x00ff0000) >> 8) |                        \
   (((x)&0x0000ff00) << 8) | (((x)&0x000000ff) << 24))

#define bswap_16(x) ((((x)&0xff00) >> 8) | (((x)&0x00ff) << 8))

enum ENDIANNESS { LT_ENDIAN, BG_ENDIAN };

inline ENDIANNESS GetEndian() {
  long int endian = 0x0000000000000001;
  return (*(char *)&endian == 0x01) ? LT_ENDIAN : BG_ENDIAN;
}

inline uint64_t htof64(uint64_t host_integer) {
  if (GetEndian() == LT_ENDIAN) {
    // Same endianness, so just return the same integer
    return host_integer;
  } else {
    // need to reverse here
    return bswap_64(host_integer);
  }
}

inline uint64_t ftoh64(uint64_t file_integer) {
  if (GetEndian() == LT_ENDIAN) {
    // Same endianness, so just return the same integer
    return file_integer;
  } else {
    // need to reverse order here
    return bswap_64(file_integer);
  }
}

inline uint32_t htof32(uint32_t host_integer) {
  if (GetEndian() == LT_ENDIAN) {
    return host_integer;
  } else {
    return bswap_32(host_integer);
  }
}

inline uint32_t ftoh32(uint32_t file_integer) {
  if (GetEndian() == LT_ENDIAN) {
    return file_integer;
  } else {
    return bswap_32(file_integer);
  }
}

inline uint16_t htof16(uint32_t host_integer) {
  if (GetEndian() == LT_ENDIAN) {
    return host_integer;
  } else {
    return bswap_16(host_integer);
  }
}

inline uint16_t ftoh16(uint16_t file_integer) {
  if (GetEndian() == LT_ENDIAN) {
    return file_integer;
  } else {
    return bswap_16(file_integer);
  }
}
