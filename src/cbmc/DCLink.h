/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINK_H
#define DCLINK_H

#include "BasicTypes.h"
#include "DCComponent.h"
#include "DCData.h"
#include "DCFactory.h"

namespace mol_setup
{
struct MolKind;
}

namespace cbmc
{
class DCLink : public DCComponent
{
public:
  virtual void PrepareNew(TrialMol& newMol, uint molIndex);
  virtual void PrepareOld(TrialMol& oldMol, uint molIndex);
  virtual void BuildOld(TrialMol& oldMol, uint molIndex);
  virtual void BuildNew(TrialMol& newMol, uint molIndex);

  virtual DCComponent* Clone()
  {
    return new DCLink(*this);
  };

  DCLink(DCData* data, const mol_setup::MolKind kind,
         uint atom, uint focus);


private:
  void IncorporateOld(TrialMol& oldMol, uint molIndex);
  void IncorporateNew(TrialMol& newMol, uint molIndex);
  void AlignBasis(TrialMol& mol);
  double GenerateDihedralsNew(TrialMol& newMol, uint molIndex,
                              double* angles, double* angleEnergy,
                              double* angleWeights);
  double GenerateDihedralsOld(TrialMol& oldMol, uint molIndex,
                              double* angles, double* angleEnergy,
                              double* angleWeights);
  void UseOldDih(TrialMol& oldMol, uint molIndex, double& energy,
                 double& weight);
  void SetOldMolBond(const uint i, const double distSq);

  DCData* data;

  bool angleFix;

  uint atom, focus, prev, prevprev;
  uint bondKind, angleKind, dihKind;
  double bondLength;
  //used for calculating 1-3 and 1-4 distance
  double bond[3];
  double oldBond[3];

  double theta, phi, thetaFix;

  double oldBondEnergy, bendEnergy, bendWeight, oneThree;
};

}

#endif
