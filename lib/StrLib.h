/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef STR_LIB_H
#define STR_LIB_H

#include <cctype>
#include <iostream> //for find_first...
#include <iterator> //For adv. version of word counting
#include <string>   //The obj. this utility collection primarily deals with
#include <vector>   //For getting a non-constant string from const vector.

namespace str {

// Converts an STL string to uppercase.
// See:
// http://stackoverflow.com/questions/735204/
// convert-a-string-in-c-to-upper-case
inline std::string ToUpper(std::string const &s) {
  std::string ret(s.size(), char());
  for (unsigned int i = 0; i < s.size(); ++i)
    ret[i] = (s[i] <= 'z' && s[i] >= 'a') ? s[i] - ('a' - 'A') : s[i];
  return ret;
}

// Wrapper/shorthand for string comparison, with support for
// insenstive comparisons.
inline bool compare(std::string const &s1, std::string const &s2,
                    const bool insensitive = true) {
  return (insensitive ? ToUpper(s1).compare(ToUpper(s2)) : s1.compare(s2)) == 0;
}

// Check if a string only holds whitespace characters.
inline bool AllWS(std::string const &str) {
  bool wht = true;
  for (unsigned int s = 0; s < str.length() && wht; s++)
    wht &= !(std::isspace(str[s]) == 0);
  return wht;
}

inline char *GetNonConstStr(std::string const &str) {
  // return non-constant copy of s.c_str()
  static std::vector<char> var;
  var.assign(str.begin(), str.end());
  var.push_back('\0');
  return &var[0];
}

inline std::string
TrimCharArrNoNullTerm(char const *const cStr, const unsigned int len,
                      const std::string &whitespace = " \t") {
  std::string str(len, ' ');
  for (unsigned int c = 0; c < len; ++c)
    str[c] = cStr[c];
  size_t strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content

  size_t strEnd = str.find_last_not_of(whitespace);
  size_t strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

inline std::string MoveTypetoStr(uint moveType) {
  std::string moveTypeStr;
  switch (moveType) {
  case 0:
    moveTypeStr = "Displacement";
    break;
  case 1:
    moveTypeStr = "Rotation";
    break;
  case 2:
    moveTypeStr = "MultiParticle";
    break;
  case 3:
    moveTypeStr = "Brownian-like MultiParticle";
    break;
  case 4:
    moveTypeStr = "Intra Molecule Transfer";
    break;
  case 5:
    moveTypeStr = "Regrowth";
    break;
  case 6:
    moveTypeStr = "Intra Molecule Exchange";
    break;
  case 7:
    moveTypeStr = "Crankshaft";
    break;
  case 8:
    moveTypeStr = "Intra Targeted Transfer";
    break;
#if ENSEMBLE == NPT
  case 9:
    moveTypeStr = "Volume Transfer";
    break;
#elif ENSEMBLE == GCMC || ENSEMBLE == GEMC
  case 9:
    moveTypeStr = "Molecule Exchange";
    break;
  case 10:
    moveTypeStr = "Molecule Transfer";
    break;
  case 11:
    moveTypeStr = "Nonequilibrium Molecule Transfer";
    break;
  case 12:
    moveTypeStr = "Targeted Transfer";
    break;
#if ENSEMBLE == GEMC
  case 13:
    moveTypeStr = "Volume Transfer";
    break;
#endif
#endif
  default:
    moveTypeStr =
        "Update MoveTypetoStr() function in lib/StrLib.h!!! Undefined";
  }

  return moveTypeStr;
}

} // namespace str

#endif /*STR_LIB_H*/
