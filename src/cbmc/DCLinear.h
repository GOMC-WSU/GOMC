/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINEAR_H
#define DCLINEAR_H
#include "CBMC.h"
#include "DCData.h"
#include <vector>

class System;
class Forcefield;
class MoleculeKind;
class Setup;

namespace cbmc
{
class DCComponent;

class DCLinear : public CBMC
{
public:
  DCLinear(System& sys, const Forcefield& ff,
           const MoleculeKind& kind, const Setup& set);

  void Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
  void Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
  void CrankShaft(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
  void BuildIDNew(TrialMol& newMol, uint molIndex);
  void BuildIDOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildGrowNew(TrialMol& newMol, uint molIndex);
  void BuildGrowOld(TrialMol& oldMol, uint molIndex);
  ~DCLinear();

private:
  uint atomSize;
  //used for when number of atom < 3
  std::vector<DCComponent*> forward, backward;
  DCComponent* idExchange;
  DCData data;
};
}

#endif
