/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCFREEHEDRONSEED_H
#define DCFREEHEDRONSEED_H
#include "DCComponent.h"
#include "DCHedron.h"
#include "CBMC.h"

namespace mol_setup { class MolKind; }

namespace cbmc {
   class DCData;

   class DCFreeHedronSeed : public DCComponent
   {  
   public:
      DCFreeHedronSeed(DCData* data, const mol_setup::MolKind& kind,
		       uint focus, uint prev);
      void PrepareNew(TrialMol& newMol, uint molIndex);
      void PrepareOld(TrialMol& oldMol, uint molIndex);
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCFreeHedronSeed(*this); };

   private:
      DCData* data;
      DCHedron hed;
      double anchorBond;
      //anchor bond energy of old molecule
      double oldBondEnergy;
      uint anchorBondKind;
   };
}

#endif
