/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
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
  void ChooseTorsion(TrialMol& mol, uint molIndex, double prevPhi[],
                     RotationMatrix& cross, RotationMatrix& tensor);
  double EvalLJ(TrialMol& mol, uint molIndex);
  void SetBasis(TrialMol& mol, uint p1, uint p2);
  void OldThetaAndPhi(TrialMol& mol, const uint atom, const uint lastAtom,
                      double& theta, double& phi) const;
  DCData* data;
  DCHedronCycle hed;
  uint nPrevBonds;
  uint prevBonded[MAX_BONDS];
  //kind[bonded][previous]
  uint dihKinds[MAX_BONDS][MAX_BONDS];

  //bond energy of built branch
  double bondEnergy;
  //bond length of prev bonded to focus
  double anchorBond, anchorBondOld;
  //bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];
  //bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
  //To find the theta and phi in bCoords
  RotationMatrix growthToWorld;
  RotationMatrix worldToGrowth;
  XYZ basisPoint;
};
}
#endif
