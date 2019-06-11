/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCFREEHEDRON_H
#define DCFREEHEDRON_H
#include "DCComponent.h"
#include "DCSingle.h"
#include "DCHedron.h"
#include "CBMC.h"

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;

class DCFreeHedron : public DCComponent
{
public:
  DCFreeHedron(DCData* data, const mol_setup::MolKind& kind,
               uint focus, uint prev);
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  void SetBondLengthNew(TrialMol& newMol);
  void SetBondLengthOld(TrialMol& oldMol);

  DCComponent* Clone()
  {
    return new DCFreeHedron(*this);
  };

private:
  DCData* data;
  DCSingle seed;
  DCHedron hed;
  //bond length of prev bonded to focus
  real anchorBond, anchorBondOld;
  //bond energy of built branch
  real bondEnergy;
  uint anchorKind;

  //bond length of atom bonded to focus
  real bondLength[MAX_BONDS];
  real bondLengthOld[MAX_BONDS];
  //bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];
};
}

#endif
