/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCLINKEDHEDRON_H
#define DCLINKEDHEDRON_H
#include "CBMC.h"
#include "DCComponent.h"
#include "DCHedron.h"
#include "TransformMatrix.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;
class DCLinkedHedron : public DCComponent {
public:
  DCLinkedHedron(DCData *data, const mol_setup::MolKind &kind, uint focus,
                 uint prev);
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void SetBondLengthNew(TrialMol &newMol);
  void SetBondLengthOld(TrialMol &oldMol);

  DCComponent *Clone() { return new DCLinkedHedron(*this); };

private:
  void ChooseTorsion(TrialMol &mol, uint molIndex, double prevPhi[],
                     RotationMatrix &cross, RotationMatrix &tensor);
  double EvalLJ(TrialMol &mol, uint molIndex);
  DCData *data;
  DCHedron hed;
  uint nPrevBonds;
  uint prevBonded[MAX_BONDS];
  // kind[bonded][previous]
  uint dihKinds[MAX_BONDS][MAX_BONDS];

  // bond energy of built branch
  double bondEnergy;
  // bond length of prev bonded to focus
  double anchorBond, anchorBondOld;
  // bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];
  // bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
};
} // namespace cbmc
#endif /*DCLINKEDHEDRON_H*/
