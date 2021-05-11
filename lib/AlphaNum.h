/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#include <vector>
#include <algorithm>
#include <cassert>
#include <sstream>      // std::stringstream
#include <algorithm>
#include <string>
#include <cctype>

struct icompare_char {
  bool operator()(char c1, char c2) {
    return std::toupper(c1) < std::toupper(c2);
  }
};

struct compare {
  bool operator()(const std::pair<std::string, int>& lhs, const std::pair<std::string, int>& rhs) {
    if (lhs.first.length() > rhs.first.length())
      return false;
    if (lhs.first.length() < rhs.first.length())
      return true;
    return std::lexicographical_compare(lhs.first.begin(), lhs.first.end(),
                                        rhs.first.begin(), rhs.first.end(),
                                        icompare_char());
  }
};

class AlphaNum
{
    public:
        AlphaNum();
        std::string uint2String(uint stringSuffix);
        uint string2Uint(std::string stringSuffix);
        struct icompare_char;
        struct compare;
};