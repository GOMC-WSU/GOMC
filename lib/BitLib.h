/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef BIT_LIB_H
#define BIT_LIB_H

#include "BasicTypes.h"
#include <vector> //for mask list

namespace bits {
static inline uint Check(const uint v, const uint pos) {
  return (v & (1 << (pos)));
}
static inline uint CountSet(const uint v) {
  uint count = 0;
  for (uint i = 0; i < 32; i++)
    if (Check(v, i))
      count++;
  return count;
}

inline std::vector<std::vector<uint>> GetMasks(uint N) {
  std::vector<std::vector<uint>> mask;
  mask.resize(N);
  for (uint i = 0; i < (uint)(1 << N) - 1; i++)
    mask[CountSet(i)].push_back(i);
  return mask;
}
} // namespace bits

#endif /*BIT_LIB_H*/
