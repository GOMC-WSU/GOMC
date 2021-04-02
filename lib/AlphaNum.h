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

class AlphaNum
{
    public:
        AlphaNum();
        std::string uint2String(uint stringSuffix);
};