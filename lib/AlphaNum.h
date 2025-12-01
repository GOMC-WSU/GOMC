/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef ALPHA_NUM_H
#define ALPHA_NUM_H

#include "BasicTypes.h" // uint
#include <algorithm>
#include <cassert>
#include <cctype>
#include <sstream> // std::stringstream
#include <string>
#include <vector>

struct icompare_char {
  bool operator()(char c1, char c2) {
    return std::toupper(c1) < std::toupper(c2);
  }
};

struct compare {
  bool operator()(const std::pair<std::string, int> &lhs,
                  const std::pair<std::string, int> &rhs) {
    if (lhs.first.length() > rhs.first.length())
      return false;
    if (lhs.first.length() < rhs.first.length())
      return true;
    return std::lexicographical_compare(lhs.first.begin(), lhs.first.end(),
                                        rhs.first.begin(), rhs.first.end(),
                                        icompare_char());
  }
};

class AlphaNum {
public:
  AlphaNum();
  std::string uint2String(uint stringSuffix);
  uint string2Uint(std::string stringSuffix);
  struct icompare_char;
  struct compare;
};

#endif /*ALPHA_NUM_H*/
