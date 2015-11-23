#ifndef DCONSPHERE_H
#define DCONSPHERE_H
#include "DCComponent.h"
#include "../../lib/BasicTypes.h"

namespace mol_setup { class MolKind; }

namespace cbmc
{
   class DCData;

   class DCOnSphere : public DCComponent
   {
    public:
      DCOnSphere(DCData* data, const mol_setup::MolKind kind,
		 uint atom, uint focus);
      void PrepareNew() {};
      void PrepareOld() {};
      void BuildOld(TrialMol& oldMol, uint molIndex);
      void BuildNew(TrialMol& newMol, uint molIndex);
      DCComponent* Clone() { return new DCOnSphere(*this); };
   private:
      DCData* data;
      uint atom, focus;
      double bondLength;
   };



}
#endif
