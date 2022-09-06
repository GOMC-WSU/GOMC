/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCCRANKSHAFTANG_H
#define DCCRANKSHAFTANG_H

#include <vector>

#include "BasicTypes.h"
#include "DCComponent.h"
#include "DCData.h"
#include "MolSetup.h"
#include "TransformMatrix.h"

using namespace mol_setup;

namespace cbmc {
class DCCrankShaftAng;

class DCCrankShaftAng : public DCComponent {
public:
  DCCrankShaftAng(DCData *data, const mol_setup::MolKind &kind, uint a0,
                  uint a1, uint a2);
  ~DCCrankShaftAng() { delete[] multiPosRotions; }
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  DCComponent *Clone() { return new DCCrankShaftAng(*this); }

private:
  void ChooseTorsion(TrialMol &mol, uint molIndex, RotationMatrix &cross,
                     RotationMatrix &tensor);
  void ChooseTorsionOld(TrialMol &mol, uint molIndex, RotationMatrix &cross,
                        RotationMatrix &tensor);
  double CalcIntraBonded(TrialMol &mol, uint molIndex);
  void ParticleNonbonded(cbmc::TrialMol const &mol, XYZArray const &trialPos,
                         const uint partIndex, const uint trials);
  void ParticleNonbonded1_N(cbmc::TrialMol const &mol, XYZArray const &trialPos,
                            const uint partIndex, const uint trials);
  void ParticleNonbonded1_4(cbmc::TrialMol const &mol, XYZArray const &trialPos,
                            const uint partIndex, const uint trials);
  void ParticleNonbonded1_3(cbmc::TrialMol const &mol, XYZArray const &trialPos,
                            const uint partIndex, const uint trials);

  DCData *data;
  XYZArray *multiPosRotions;
  uint a0, a1, a2, numAtom, totAtoms;
  std::vector<uint> atoms;
  std::vector<Angle> ang;
  std::vector<Dihedral> dih;
};
} // namespace cbmc
#endif
