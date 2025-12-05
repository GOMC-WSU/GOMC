/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCCOMPONENT_H
#define DCCOMPONENT_H
#include "BasicTypes.h"

// Class for Deferred Coupling CBMC components

namespace cbmc {
class TrialMol;

class DCComponent {
public:
  // Perform Decoupled portions of CBMC
  virtual void PrepareNew(TrialMol &newMol, uint molIndex) {}
  virtual void PrepareOld(TrialMol &oldMol, uint molIndex) {}

  // Perform Coupled final build
  virtual void BuildOld(TrialMol &oldMol, uint molIndex) = 0;
  virtual void BuildNew(TrialMol &newMol, uint molIndex) = 0;

  virtual void UpdateAcceptance(const TrialMol &mol) {}
  virtual ~DCComponent(){};
};
} // namespace cbmc

#endif /*DCCOMPONENT_H*/
