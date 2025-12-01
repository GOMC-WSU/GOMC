/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCSINGLE_H
#define DCSINGLE_H

#include "BasicTypes.h"
#include "DCComponent.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;

class DCSingle : public DCComponent {
public:
  DCSingle(DCData *data, uint atom);
  void PrepareNew(TrialMol &newMol, uint molIndex) {};
  void PrepareOld(TrialMol &oldMol, uint molIndex) {};
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  DCComponent *Clone() { return new DCSingle(*this); }

private:
  DCData *data;
  uint atom;
};
} // namespace cbmc
#endif /*DCSINGLE_H*/
