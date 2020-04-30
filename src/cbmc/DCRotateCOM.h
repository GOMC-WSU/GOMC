/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCROTATECOM_H
#define DCROTATECOM_H
#include "DCComponent.h"
#include "XYZArray.h"
#include "CBMC.h"
#include "GeomLib.h"

namespace mol_setup
{
class MolKind;
}

namespace cbmc
{
class DCData;

class DCRotateCOM : public DCComponent
{
public:
  DCRotateCOM(DCData* data, const mol_setup::MolKind kind);
  ~DCRotateCOM()
  {
    delete[] multiPosRotions;
  }
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void PickTransferCOMNew(TrialMol& newMol, uint molIndex);
  void PickTransferCOMOld(TrialMol& oldMol, uint molIndex);
  void BuildOld(TrialMol& oldMol, uint molIndex);
  void BuildNew(TrialMol& newMol, uint molIndex);
  void RandRotateZ();
  DCComponent* Clone()
  {
    return new DCRotateCOM(*this);
  };

private:
  DCData* data;
  XYZ COM, oldCOM;
  //rotation matrix around z axis
  XYZArray rotateMatrix;
  //inverse of matrix
  XYZArray invMatrix;
  XYZArray *multiPosRotions;
  uint atomNumber;
};
}

#endif
