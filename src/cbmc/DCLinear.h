/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCLINEAR_H
#define DCLINEAR_H
#include <vector>

#include "CBMC.h"
#include "DCData.h"

class System;
class Forcefield;
class MoleculeKind;
class Setup;

namespace cbmc {
class DCComponent;

class DCLinear : public CBMC {
public:
  DCLinear(System &sys, const Forcefield &ff, const MoleculeKind &kind,
           const Setup &set);

  void Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void BuildIDNew(TrialMol &newMol, uint molIndex);
  void BuildIDOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildGrowNew(TrialMol &newMol, uint molIndex);
  void BuildGrowOld(TrialMol &oldMol, uint molIndex);
  // used in TargetedSwap
  void BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  ~DCLinear();

private:
  uint atomSize;
  // used for when number of atom < 3
  std::vector<DCComponent *> forward, backward;
  DCComponent *idExchange;
  DCData data;
};
} // namespace cbmc

#endif /*DCLINEAR_H*/
