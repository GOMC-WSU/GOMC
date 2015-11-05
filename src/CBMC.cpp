/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "CBMC.h"
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "cbmc/DCLinear.h"
#include "cbmc/DCGraph.h"
#include <vector>


namespace cbmc{

   CBMC* MakeCBMC(System& sys, const Forcefield& ff,
                  const MoleculeKind& kind, const Setup& set)
   {
      std::vector<uint> bondCount(kind.NumAtoms(), 0);
      for (uint i = 0; i < kind.bondList.count; ++i)
      {
         bondCount[kind.bondList.part1[i]]++;
         bondCount[kind.bondList.part2[i]]++;
      }
      bool branched = false;
      for (uint i = 0; i < kind.NumAtoms(); ++i)
      {
         if (bondCount[i] > 2)
            branched = true;
      }
      bool cyclic = (kind.NumBonds() > kind.NumAtoms() - 1) ? true : false;

      if (branched)
      {
         return new DCGraph(sys, ff, kind, set);
      }
      else {
         return new DCLinear(sys, ff, kind, set);
      }
   }

}
