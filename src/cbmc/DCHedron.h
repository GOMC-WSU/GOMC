/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCHEDRON_H
#define DCHEDRON_H
#include "DCComponent.h"
#include "CBMC.h"

namespace mol_setup { class MolKind; }

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
      uint Bonded(uint i) const { return bonded[i]; }
      double Theta(uint i) const { return theta[i]; }
      double Phi(uint i) const { return phi[i]; }
      double BondLength(uint i) const { return newBondLength[i]; }
      double BondLengthOld(uint i) const { return oldBondLength[i]; }
      double GetNewAnchor() const {return newAnchorBond; }
      double GetOldAnchor() const {return oldAnchorBond; }
      double GetWeight();
      double GetEnergy() { return bendEnergy; }
      double GetNonBondedEn() { return oneThree; }
      double GetOldBondEn() { return oldBondEnergy; }
      double GetNewBondEn() { return newBondEnergy; }
      double GetOldBondW() { return oldBondWeight; }
      double GetNewBondW() { return newBondWeight; }
      uint NumBond() { return nBonds; }
      uint Focus() { return focus; }
      uint Prev() { return prev; }
      
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
      void SetOldBond(TrialMol& oldMol); 
      void SetNewBond(TrialMol& newMol);

      DCData* data;
      uint focus, prev;
      uint nBonds;
      //atoms bonded to focus, being build

      //bond length of atom bonded to focus
      double eqBondLength[MAX_BONDS];
      double newBondLength[MAX_BONDS];
      double oldBondLength[MAX_BONDS];
      
      //angleKinds[i][j] = kind between bonded[i] and bonded[j]
      //except angleKinds[i][i] = kind between bonded[i] and prev
      uint angleKinds[MAX_BONDS][MAX_BONDS];
      //bondKind between bonded[i] and focus
      uint bondKinds[MAX_BONDS];
      uint anchorKind;

      double theta[MAX_BONDS];
      double thetaWeight[MAX_BONDS];
      double phi[MAX_BONDS];
      double phiWeight[MAX_BONDS];
      double bendEnergy, oneThree;
      double eqAnchorBond, newAnchorBond, oldAnchorBond;
      double oldBondEnergy, oldBondWeight;
      double newBondEnergy, newBondWeight;
   };
}

#endif

