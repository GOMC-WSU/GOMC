/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCONSPHERE_H
#define DCONSPHERE_H
#include "BasicTypes.h"
#include "DCComponent.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;

class DCOnSphere : public DCComponent {
public:
  DCOnSphere(DCData *data, const mol_setup::MolKind kind, uint atom,
             uint focus);
  void PrepareNew(TrialMol &newMol, uint molIndex) {};
  void PrepareOld(TrialMol &oldMol, uint molIndex) {};
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void SetBondLengthNew(TrialMol &newMol);
  void SetBondLengthOld(TrialMol &oldMol);

  DCComponent *Clone() { return new DCOnSphere(*this); };

private:
  double BondEnergyNew(TrialMol &newMol);
  double BondEnergyOld(TrialMol &oldMol);
  DCData *data;
  uint atom, focus;
  uint bondKind;
  double bondLength, bondLengthOld;
  double bondEnergy;
};

} // namespace cbmc
#endif /*DCONSPHERE_H*/
