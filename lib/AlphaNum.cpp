/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "AlphaNum.h"

AlphaNum::AlphaNum() {}

/* Alphanumeric A, B, ... Z, AA, AB, ... AZ, */
/* Index from 0 to A */
// return true if s1 comes before s2

// compare character case-insensitive

std::string AlphaNum::uint2String(uint stringSuffix) {

  std::stringstream ss;
  int intermediate = stringSuffix;
  do {
    char charSuffix = 'A';
    int remainder = intermediate % 26;
    /* Increment Char A until reach suffix or 27 which will be Z. */
    for (int j = 0; j < remainder; j++) {
      charSuffix++;
    }
    ss << charSuffix;
    intermediate /= 26;
    intermediate--;
  } while (intermediate >= 0);

  std::string backwards = ss.str();
  std::string forwards;

  return forwards.assign(backwards.rbegin(), backwards.rend());
}

uint AlphaNum::string2Uint(const std::string &stringSuffix) {

  const char charSuffix = 'A';
  uint index = 0;
  for (int i = 0; i < stringSuffix.length(); ++i) {
    index += 26 * i + (stringSuffix[i] - charSuffix);
  }
  return index;
}
