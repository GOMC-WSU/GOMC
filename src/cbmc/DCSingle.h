/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCSINGLE_H
#define DCSINGLE_H

#include "DCComponent.h"
#include "BasicTypes.h"

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;

class DCSingle : public DCComponent
{
public:
  DCSingle(DCData* data, uint atom) : data(data), atom(atom) {}
  void PrepareNew(TrialMol& newMol, uint molIndex) {};
  void PrepareOld(TrialMol& oldMol, uint molIndex) {};
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  DCComponent* Clone()
  {
    return new DCSingle(*this);
  }

private:
  DCData* data;
  uint atom;
};
}
#endif
