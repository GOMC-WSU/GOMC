#ifndef DCLINKNODIH_H
#define DCLINKNODIH_H

#include "../../lib/BasicTypes.h"
#include "DCComponent.h"
#include "DCData.h"

namespace cbmc
{
   class DCLinkNoDih : public DCComponent
   {
    public:
      DCLinkNoDih(DCData* data, const mol_setup::MolKind kind,
		  uint atom, uint focus);

      virtual void PrepareNew();
      virtual void PrepareOld();
      virtual void BuildOld(TrialMol& oldMol, uint molIndex);
      virtual void BuildNew(TrialMol& newMol, uint molIndex);

      virtual DCLinkNoDih* Clone() { return new DCLinkNoDih(*this); };

    private:
      void IncorporateOld(TrialMol& oldMol);
      void IncorporateNew(TrialMol& newMol);
      void AlignBasis(TrialMol& mol);

      DCData* data;
      uint atom, focus, prev;
      uint angleKind;
      double bondLength;
      double theta;
      double bendEnergy, bendWeight;
   };
}

#endif
