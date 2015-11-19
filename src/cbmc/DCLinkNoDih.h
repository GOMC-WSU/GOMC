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

      virtual void PrepareNew(TrialMol& newMol, uint molIndex);
      virtual void PrepareOld(TrialMol& oldMol, uint molIndex);
      virtual void BuildOld(TrialMol& oldMol, uint molIndex);
      virtual void BuildNew(TrialMol& newMol, uint molIndex);

      virtual DCLinkNoDih* Clone() { return new DCLinkNoDih(*this); };

    private:
      void IncorporateOld(TrialMol& oldMol, uint molIndex);
      void IncorporateNew(TrialMol& newMol, uint molIndex);
      void AlignBasis(TrialMol& mol);
      void SetOldMolBond(const uint i, const double distSq);

      DCData* data;
      uint atom, focus, prev;
      uint bondKind, angleKind;
      double bondLength;
      double bond[2];
      double oldBond[2];
      double theta;
      double oldBondEnergy, bendEnergy, bendWeight, oneThree;
   };
}

#endif
