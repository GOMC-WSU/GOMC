/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BIT_LIB_H
#define BIT_LIB_H

#include <vector> //for mask list
#include "BasicTypes.h"

namespace bits
{
   static inline uint Check(const uint v, const uint pos) 
   { return (v & (1<<(pos))); }
   static inline uint CountSet(const uint v)
   { 
      uint count=0;
      for (uint i = 0; i < 32; i++)
	 if (Check(v,i)) 
	   count++;
      return count;
   }

   inline std::vector< std::vector<uint> > GetMasks(uint N)
   {
      std::vector< std::vector<uint> > mask;
      mask.resize(N);
      for (uint i = 0; i < (uint)(1<<N)-1; i++)
	 mask[CountSet(i)].push_back(i);
      return mask;
   }
}

#endif /*BIT_LIB_H*/

