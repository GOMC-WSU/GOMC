/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCFREEHEDRONSEED_H
#define DCFREEHEDRONSEED_H
#include "CBMC.h"
#include "DCComponent.h"
#include "DCHedron.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;

class DCFreeHedronSeed : public DCComponent {
public:
  DCFreeHedronSeed(DCData *data, const mol_setup::MolKind &kind, uint focus,
                   uint prev);
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void SetBondLengthNew(TrialMol &newMol);
  void SetBondLengthOld(TrialMol &oldMol);

  DCComponent *Clone() { return new DCFreeHedronSeed(*this); };

private:
  DCData *data;
  DCHedron hed;
  // bond length of prev bonded to focus
  double anchorBond, anchorBondOld;
  // bond energy of built branch
  double bondEnergy;
  uint anchorKind;

  // bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];
  // bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
};
} // namespace cbmc

#endif
