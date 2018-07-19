/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCCRANKSHAFTANG_H
#define DCCRANKSHAFTANG_H

#include "DCComponent.h"
#include "DCData.h"
#include "BasicTypes.h"
#include "TransformMatrix.h"
#include <vector>

namespace mol_setup
{
  class MolKind;
}

namespace cbmc
{
class DCCrankShaftAng;

class DCCrankShaftAng : public DCComponent
{
public:
  DCCrankShaftAng(DCData* data, const mol_setup::MolKind& kind,
                  uint a0, uint a1, uint a2);
  ~DCCrankShaftAng() {delete[] multiPosRotions;}
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  DCComponent* Clone()
  {
    return new DCCrankShaftAng(*this);
  }

private:
  DCData* data;
  XYZArray *multiPosRotions;
  uint a0, a1, a2, numAtom, totAtoms;
  std::vector<uint> atoms;
};
}
#endif
