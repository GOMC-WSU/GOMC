#ifndef DCHEDRON_H
#define DCHEDRON_H
#include "DCComponent.h"
#include "../CBMC.h"

namespace mol_setup
{
struct MolKind;
}

namespace cbmc
{
class DCData;

class DCHedron
{
public:
  DCHedron(DCData* data, const mol_setup::MolKind& kind,
           uint focus, uint prev);
  void PrepareNew(TrialMol& newMol, uint molIndex);
  void PrepareOld(TrialMol& oldMol, uint molIndex);
  void IncorporateOld(TrialMol& oldMol, uint molIndex);
  void ConstrainedAnglesOld(uint nTrials, TrialMol& oldMol, uint molIndex);
  uint Bonded(uint i) const
  {
    return bonded[i];
  }
  double Theta(uint i) const
  {
    return theta[i];
  }
  double Phi(uint i) const
  {
    return phi[i];
  }
  double BondLength(uint i) const
  {
    return bondLength[i];
  }
  double BondLengthOld(uint i) const
  {
    return bondLengthOld[i];
  }
  double GetWeight();
  double GetEnergy()
  {
    return bendEnergy;
  }
  double GetNonBondedEn()
  {
    return oneThree;
  }
  double GetOldBondEn()
  {
    return oldBondEnergy;
  }
  uint NumBond()
  {
    return nBonds;
  }
  uint Focus()
  {
    return focus;
  }
  uint Prev()
  {
    return prev;
  }

  //need to go to private
  uint bonded[MAX_BONDS];


private:
  void GenerateAnglesNew(TrialMol& newMol, uint molIndex, uint kind,
                         uint nTrials, uint bType);
  void GenerateAnglesOld(TrialMol& oldMol, uint molIndex, uint kind,
                         uint nTrials, uint bType);
  void FreeAnglesNew(TrialMol& newMol, uint molIndex, uint nTrials);
  void FreeAnglesOld(TrialMol& oldMol, uint molIndex, uint nTrials);
  void ConstrainedAngles(TrialMol& newMol, uint molIndex, uint nTrials);


  DCData* data;
  uint focus, prev;
  uint nBonds;
  //atoms bonded to focus, being build

  //bond length of atom bonded to focus
  double bondLength[MAX_BONDS];
  double bondLengthOld[MAX_BONDS];

  //angleKinds[i][j] = kind between bonded[i] and bonded[j]
  //except angleKinds[i][i] = kind between bonded[i] and prev
  uint angleKinds[MAX_BONDS][MAX_BONDS];
  //bondKind between bonded[i] and focus
  uint bondKinds[MAX_BONDS];

  double theta[MAX_BONDS];
  double thetaWeight[MAX_BONDS];
  double phi[MAX_BONDS];
  double phiWeight[MAX_BONDS];
  double bendEnergy, oldBondEnergy, oneThree;
  double anchorBond, anchorBondOld;
};
}

#endif

