/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCROTATECOM_H
#define DCROTATECOM_H
#include "CBMC.h"
#include "DCComponent.h"
#include "GeomLib.h"
#include "XYZArray.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;

class DCRotateCOM : public DCComponent {
public:
  DCRotateCOM(DCData *data, const mol_setup::MolKind kind);
  ~DCRotateCOM() { delete[] multiPosRotions; }
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void PickTransferCOMNew(TrialMol &newMol, uint molIndex);
  void PickTransferCOMOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void RandRotateZ();
  DCComponent *Clone() { return new DCRotateCOM(*this); };

private:
  DCData *data;
  XYZ COM, oldCOM;
  // rotation matrix around z axis
  XYZArray rotateMatrix;
  // inverse of matrix
  XYZArray invMatrix;
  XYZArray *multiPosRotions;
  uint atomNumber;
};
} // namespace cbmc

#endif /*DCROTATECOM_H*/
