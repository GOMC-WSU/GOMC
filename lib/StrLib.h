/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef STR_LIB_H
#define STR_LIB_H

#include <string> //The obj. this utility collection primarily deals with
#include <vector> //For getting a non-constant string from const vector.
#include <iterator> //For adv. version of word counting
#include <cctype>
#include <iostream> //for find_first...

namespace str
{

// Converts an STL string to uppercase.
// See:
//http://stackoverflow.com/questions/735204/
//convert-a-string-in-c-to-upper-case
inline std::string ToUpper(std::string const& s)
{
  std::string ret(s.size(), char());
  for(unsigned int i = 0; i < s.size(); ++i)
    ret[i] = (s[i] <= 'z' && s[i] >= 'a') ? s[i] - ('a' - 'A') : s[i];
  return ret;
}

// Wrapper/shorthand for string comparison, with support for
// insenstive comparisons.
inline bool compare(std::string const& s1, std::string const& s2,
                    const bool insensitive = true)
{
  return (insensitive ? ToUpper(s1).compare(ToUpper(s2)) : s1.compare(s2)) == 0;
}

//Check if a string only holds whitespace characters.
inline bool AllWS(std::string const& str)
{
  bool wht = true;
  for (unsigned int s = 0; s < str.length() && wht; s++)
    wht &= !(std::isspace(str[s]) == 0);
  return wht;
}

inline char* GetNonConstStr(std::string const& str)
{
  //return non-constant copy of s.c_str()
  static std::vector<char> var;
  var.assign(str.begin(), str.end());
  var.push_back('\0');
  return &var[0];
}

inline std::string TrimCharArrNoNullTerm
(char const*const cStr,
 const unsigned int len,
 const std::string& whitespace = " \t")
{
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

inline std::string MoveTypetoStr(uint moveType)
{
  std::string moveTypeStr;
  switch (moveType) {
    case 0:
      moveTypeStr = "DISPLACE";
      break;
    case 1:
      moveTypeStr = "ROTATE";
      break;
    case 2:
      moveTypeStr = "MULTIPARTICLE";
      break;
    case 3:
      moveTypeStr = "BROWNIAN MULTIPARTICLE";
      break;
    case 4:
      moveTypeStr = "INTRA_SWAP";
      break;
    case 5:
      moveTypeStr = "REGROWTH";
      break;
    case 6:
      moveTypeStr = "INTRA_MEMC";
      break;
    case 7:
      moveTypeStr = "CRANKSHAFT";
      break;
    case 8:
      moveTypeStr = "INTRA_TARGETED_SWAP";
      break;
#if ENSEMBLE == NPT
    case 9:
      moveTypeStr = "VOL_TRANSFER";
      break;
#elif ENSEMBLE == GCMC || ENSEMBLE == GEMC
    case 9:
      moveTypeStr = "MEMC";
      break;
    case 10:
      moveTypeStr = "MOL_TRANSFER";
      break;
    case 11:
      moveTypeStr = "NE_MTMC";
      break;
    case 12:
      moveTypeStr = "TARGETED_SWAP";
      break;
#elif ENSEMBLE == GEMC
    case 13:
      moveTypeStr = "VOL_TRANSFER";
      break;
#endif
    default:
      moveTypeStr = "Update printMove() function!!! Undefined";
  }

  return moveTypeStr;
}

}

#endif /*STR_LIB_H*/
