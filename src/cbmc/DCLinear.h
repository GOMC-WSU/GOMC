/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINEAR_H
#define DCLINEAR_H
#include "../CBMC.h"
#include "DCData.h"
#include <vector>

class System;
class Forcefield;
class MoleculeKind;
class Setup;

namespace cbmc {
   class DCComponent;

   class DCLinear : public CBMC
   {
   public:
      DCLinear(System& sys, const Forcefield& ff,
         const MoleculeKind& kind, const Setup& set);

      void Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
      ~DCLinear();

   private:
      std::vector<DCComponent*> forward, backward;
      DCData data;
   };
}

#endif
