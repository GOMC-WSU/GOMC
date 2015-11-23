#ifndef DCCOMPONENT_H
#define DCCOMPONENT_H
#include "../../lib/BasicTypes.h"

//Class for Deferred Coupling CBMC components

namespace cbmc {
   class TrialMol;

   class DCComponent
   {
   public:
      //Perform Decoupled portions of CBMC
      virtual void PrepareNew() {}
      virtual void PrepareOld() {}

      //Perform Coupled final build
      virtual void BuildOld(TrialMol& oldMol, uint molIndex) = 0;
      virtual void BuildNew(TrialMol& newMol, uint molIndex) = 0;

      virtual void UpdateAcceptance(const TrialMol& mol) {}
      virtual ~DCComponent() {};
   };
}

#endif