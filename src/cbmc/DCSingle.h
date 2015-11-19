#ifndef DCSINGLE_H
#define DCSINGLE_H

#include "DCComponent.h"
#include "../../lib/BasicTypes.h"

namespace mol_setup { struct MolKind; }

namespace cbmc
{
   class DCData;

   class DCSingle : public DCComponent
   {
    public:
      DCSingle(DCData* data, uint atom) : data(data), atom(atom) {}
      void PrepareNew(TrialMol& newMol, uint molIndex) {};
      void PrepareOld(TrialMol& oldMol, uint molIndex) {};
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCSingle(*this); }

   private:
      DCData* data;
      uint atom;
   };
}
#endif
