/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCONSPHERE_H
#define DCONSPHERE_H
#include "DCComponent.h"
#include "BasicTypes.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;

class DCOnSphere : public DCComponent
{
public:
  DCOnSphere(DCData* data, const mol_setup::MolKind kind,
             uint atom, uint focus);
  void PrepareNew(TrialMol& newMol, uint molIndex) {};
  void PrepareOld(TrialMol& oldMol, uint molIndex) {};
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  void SetBondLengthNew(TrialMol& newMol);
  void SetBondLengthOld(TrialMol& oldMol);

  DCComponent* Clone()
  {
    return new DCOnSphere(*this);
  };

private:
  real BondEnergyNew(TrialMol& newMol);
  real BondEnergyOld(TrialMol& oldMol);
  DCData* data;
  uint atom, focus;
  uint bondKind;
  real bondLength, bondLengthOld;
  real bondEnergy;
};



}
#endif
