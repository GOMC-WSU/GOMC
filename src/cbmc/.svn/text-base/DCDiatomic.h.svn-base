#ifndef DCDIATOMIC_H
#define DCDIATOMIC_H

#include "DCComponent.h"
#include "../../lib/BasicTypes.h"

namespace mol_setup { struct MolKind; }

namespace cbmc {
   class DCData;

   class DCDiatomic : public DCComponent
   {
   public:
      DCDiatomic(DCData* data, const mol_setup::MolKind kind, 
		 uint first, uint second);
      void PrepareNew() {};
      void PrepareOld() {};
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCDiatomic(*this); };

   private:
      DCData* data;
      uint first, second;
      double bondLength;
   };
}
#endif
