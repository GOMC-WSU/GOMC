/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINKEDCYCLE_H
#define DCLINKEDCYCLE_H
#include "DCComponent.h"
#include "CBMC.h"
#include "DCHedronCycle.h"
#include "TransformMatrix.h"

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;
class DCLinkedCycle : public DCComponent
{
public:
  DCLinkedCycle(DCData* data, const mol_setup::MolKind& kind,
                std::vector<int> cycAtoms, uint focus, uint prev);
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  void SetBondLengthNew(TrialMol& newMol);
  void SetBondLengthOld(TrialMol& oldMol);

  DCComponent* Clone()
  {
    return new DCLinkedCycle(*this);
  };

private:
  void ChooseTorsion(TrialMol& mol, uint molIndex, real prevPhi[],
                     RotationMatrix& cross, RotationMatrix& tensor);
  real EvalLJ(TrialMol& mol, uint molIndex);
  //Calculate the dihedral using bCoords
  real CalcDih(TrialMol& mol, uint a0, uint a1, uint a2, uint a3);
  void CaclIntraEnergy(TrialMol& mol, const uint bIdx, const uint molIndex);

  DCData* data;
  DCHedronCycle hed;
  uint nPrevBonds;
  uint prevBonded[MAX_BONDS];
  //kind[bonded][previous]
  uint dihKinds[MAX_BONDS][MAX_BONDS];
  //Used in finding the atom bonded to prev and focus and bith are in the ring
  uint prevBondedRing, focBondedRing;
  //Calculate torsion difference to match ring dihedral
  real torDiff;

  //bond energy of built branch
  real bondEnergy;
  //bond length of prev bonded to focus
  real anchorBond, anchorBondOld;
  //bond length of atom bonded to focus
  real bondLength[MAX_BONDS];
  real bondLengthOld[MAX_BONDS];
  //bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
  bool bondedInRing[MAX_BONDS];

  std::vector< std::vector<mol_setup::Dihedral> > bondedFocusDih;
  std::vector<bool> bExist;
};
}
#endif
