/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINKEDHEDRON_H
#define DCLINKEDHEDRON_H
#include "DCComponent.h"
#include "CBMC.h"
#include "DCHedron.h"
#include "TransformMatrix.h"

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;
class DCLinkedHedron : public DCComponent
{
public:
  DCLinkedHedron(DCData* data, const mol_setup::MolKind& kind,
                 uint focus, uint prev);
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  DCComponent* Clone()
  {
    return new DCLinkedHedron(*this);
  };
private:
  void ChooseTorsion(TrialMol& mol, uint molIndex, double prevPhi[],
		     RotationMatrix& cross, RotationMatrix& tensor);
  double EvalLJ(TrialMol& mol, uint molIndex);
  DCData* data;
  DCHedron hed;
  uint nPrevBonds;
  uint prevBonded[MAX_BONDS];
  //kind[bonded][previous]
  uint dihKinds[MAX_BONDS][MAX_BONDS];
};
}
#endif
