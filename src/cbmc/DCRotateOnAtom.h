/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCROTATEONATOM
#define DCROTATEONATOM

#include "DCComponent.h"
#include "DCData.h"
#include "BasicTypes.h"
#include "TransformMatrix.h"
#include "MolSetup.h"
#include <vector>

using namespace mol_setup;

namespace cbmc
{
class DCRotateOnAtom;

class DCRotateOnAtom : public DCComponent
{
public:
  DCRotateOnAtom(DCData* data, const mol_setup::MolKind& kind,
                 uint a0, uint a1, uint a2);
  ~DCRotateOnAtom()
  {
    delete[] multiPosRotions;
  }
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  DCComponent* Clone()
  {
    return new DCRotateOnAtom(*this);
  }

private:
  void ChooseTorsion(TrialMol& mol, uint molIndex, RotationMatrix& cross,
                     RotationMatrix& tensor);
  void ChooseTorsionOld(TrialMol& mol, uint molIndex, RotationMatrix& cross,
                        RotationMatrix& tensor);
  real CalcIntraBonded(TrialMol& mol, uint molIndex);
  void ParticleNonbonded(cbmc::TrialMol const& mol, XYZArray const& trialPos,
                         const uint partIndex, const uint trials);
  void ParticleNonbonded1_N(cbmc::TrialMol const& mol, XYZArray const& trialPos,
                            const uint partIndex, const uint trials);
  void ParticleNonbonded1_4(cbmc::TrialMol const& mol, XYZArray const& trialPos,
                            const uint partIndex, const uint trials);
  void ParticleNonbonded1_3(cbmc::TrialMol const& mol, XYZArray const& trialPos,
                            const uint partIndex, const uint trials);

  DCData* data;
  XYZArray *multiPosRotions;
  uint a0, a1, a2, numAtom, totAtoms;
  std::vector<uint> atoms;
  std::vector<Angle> ang;
  std::vector<Dihedral> dih;
};
}
#endif
