/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCLINKEDCYCLE_H
#define DCLINKEDCYCLE_H
#include "CBMC.h"
#include "DCComponent.h"
#include "DCHedronCycle.h"
#include "TransformMatrix.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;
class DCLinkedCycle : public DCComponent {
public:
  DCLinkedCycle(DCData *data, const mol_setup::MolKind &kind,
                std::vector<int> cycAtoms, uint focus, uint prev);
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void SetBondLengthNew(TrialMol &newMol);
  void SetBondLengthOld(TrialMol &oldMol);

  DCComponent *Clone() { return new DCLinkedCycle(*this); };

private:
  void ChooseTorsion(TrialMol &mol, uint molIndex, double prevPhi[],
                     RotationMatrix &cross, RotationMatrix &tensor);
  double EvalLJ(TrialMol &mol, uint molIndex);
  // Calculate the dihedral using bCoords
  double CalcDih(TrialMol &mol, uint a0, uint a1, uint a2, uint a3);
  void CaclIntraEnergy(TrialMol &mol, const uint bIdx, const uint molIndex);

  DCData *data;
  DCHedronCycle hed;
  uint nPrevBonds;
  uint prevBonded[MAX_BONDS];
  // kind[bonded][previous]
  uint dihKinds[MAX_BONDS][MAX_BONDS];
  // Used in finding the atom bonded to prev and focus and both are in the ring
  int prevBondedRing, focBondedRing;
  // Calculate torsion difference to match ring dihedral
  double torDiff;

  // bond energy of built branch
  double bondEnergy;
  // bond length of prev bonded to focus
  double anchorBond, anchorBondOld;
  // bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];
  // bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
  bool bondedInRing[MAX_BONDS];

  std::vector<std::vector<mol_setup::Dihedral>> bondedFocusDih;
  std::vector<bool> bExist;
};
} // namespace cbmc
#endif /*DCLINKEDCYCLE_H*/
