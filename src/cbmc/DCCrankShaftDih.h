/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCCRANKSHAFTDIH_H
#define DCCRANKSHAFTDIH_H

#include <vector>

#include "BasicTypes.h"
#include "DCComponent.h"
#include "DCData.h"
#include "MolSetup.h"
#include "TransformMatrix.h"

using namespace mol_setup;

namespace cbmc {
class DCCrankShaftDih;

class DCCrankShaftDih : public DCComponent {
public:
  DCCrankShaftDih(DCData *data, const mol_setup::MolKind &kind, uint a0,
                  uint a1, uint a2, uint a3);
  ~DCCrankShaftDih() { delete[] multiPosRotions; }
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  DCComponent *Clone() { return new DCCrankShaftDih(*this); }

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
  uint a0, a1, a2, a3, numAtom, totAtoms;
  std::vector<uint> atoms;
  std::vector<Angle> ang;
  std::vector<Dihedral> dih;
};
} // namespace cbmc
#endif /*DCCRANKSHAFTDIH_H*/
