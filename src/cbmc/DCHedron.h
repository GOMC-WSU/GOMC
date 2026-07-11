/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef DCHEDRON_H
#define DCHEDRON_H
#include "CBMC.h"
#include "DCComponent.h"

namespace mol_setup {
class MolKind;
}

namespace cbmc {
class DCData;

class DCHedron {
public:
  DCHedron(DCData *data, const mol_setup::MolKind &kind, uint focus, uint prev);
  void PrepareNew(TrialMol &newMol, uint molIndex);
  void PrepareOld(TrialMol &oldMol, uint molIndex);
  void IncorporateOld(TrialMol &oldMol, uint molIndex);
  void ConstrainedAnglesOld(uint nTrials, TrialMol &oldMol, uint molIndex);
  void SetBondNew(double const *bondLen, double const &anchBond);
  void SetBondOld(double const *bondLen, double const &anchBond);
  uint Bonded(uint i) const { return bonded[i]; }
  double Theta(uint i) const { return theta[i]; }
  double Phi(uint i) const { return phi[i]; }
  double GetWeight() const;
  double GetEnergy() const { return bendEnergy; }
  double GetNonBondedEn() const { return oneThree; }
  uint NumBond() const { return nBonds; }
  uint Focus() const { return focus; }
  uint Prev() const { return prev; }

  // need to go to private
  uint bonded[MAX_BONDS];

private:
  void GenerateAnglesNew(TrialMol &newMol, uint molIndex, uint kind,
                         uint nTrials, uint bType);
  void GenerateAnglesOld(TrialMol &oldMol, uint molIndex, uint kind,
                         uint nTrials, uint bType);
  void FreeAnglesNew(TrialMol &newMol, uint molIndex, uint nTrials);
  void FreeAnglesOld(TrialMol &oldMol, uint molIndex, uint nTrials);
  void ConstrainedAngles(TrialMol &newMol, uint molIndex, uint nTrials);

  DCData *data;
  uint focus, prev;
  uint nBonds;
  // atoms bonded to focus, being build

  // bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];

  // angleKinds[i][j] = kind between bonded[i] and bonded[j]
  // except angleKinds[i][i] = kind between bonded[i] and prev
  uint angleKinds[MAX_BONDS][MAX_BONDS];
  double theta[MAX_BONDS];
  double thetaWeight[MAX_BONDS];
  double phi[MAX_BONDS];
  double phiWeight[MAX_BONDS];
  double bendEnergy, oneThree;
  double anchorBond, anchorBondOld;
};
} // namespace cbmc

#endif /*DCHEDRON_H*/
