/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCHEDRON_H
#define DCHEDRON_H
#include "DCComponent.h"
#include "../CBMC.h"

namespace mol_setup { struct MolKind; }

namespace cbmc
{
   class DCData;

   class DCHedron
   {
    public:
      DCHedron(DCData* data, const mol_setup::MolKind& kind,
	       uint focus, uint prev);
      void PrepareNew();
      void PrepareOld();
      void IncorporateOld(TrialMol& oldMol);
      void ConstrainedAnglesOld(uint nTrials, TrialMol& oldMol);
      uint Bonded(uint i) const { return bonded[i]; }
      double Theta(uint i) const { return theta[i]; }
      double Phi(uint i) const { return phi[i]; }
      double BondLength(uint i) const { return bondLength[i]; }
      double GetWeight();
      double GetEnergy() { return bendEnergy; }
      uint NumBond() { return nBonds; }
      uint Focus() { return focus; }
      uint Prev() { return prev; }

    private:
      void GenerateAngles(uint kind, uint nTrials);
      void FreeAngles(uint nTrials);
      void ConstrainedAngles(uint nTrials);


      DCData* data;
      uint focus, prev;
      uint nBonds;
      //atoms bonded to focus, being build
      uint bonded[MAX_BONDS];
      //atoms bonded to prev, in dihs
      double bondLength[MAX_BONDS];
      //angleKinds[i][j] = kind between bonded[i] and bonded[j]
      //except angleKinds[i][i] = kind between bonded[i] and prev
      uint angleKinds[MAX_BONDS][MAX_BONDS];

      double theta[MAX_BONDS];
      double thetaWeight[MAX_BONDS];
      double phi[MAX_BONDS];
      double phiWeight[MAX_BONDS];
      double bendEnergy;
   };
}

#endif

