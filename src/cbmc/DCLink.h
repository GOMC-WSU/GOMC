/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCLINK_H
#define DCLINK_H

#include "../../lib/BasicTypes.h"
#include "DCComponent.h"
#include "DCData.h"
#include "DCFactory.h"

namespace mol_setup { class MolKind; }

namespace cbmc
{
   class DCLink : public DCComponent
   {
    public:
      virtual void PrepareNew();
      virtual void PrepareOld();
      virtual void BuildOld(TrialMol& oldMol, uint molIndex);
      virtual void BuildNew(TrialMol& newMol, uint molIndex);

      virtual DCComponent* Clone() { return new DCLink(*this); };

      DCLink(DCData* data, const mol_setup::MolKind kind,
	     uint atom, uint focus);


    private:
      void IncorporateOld(TrialMol& oldMol);
      void IncorporateNew(TrialMol& newMol);
      void AlignBasis(TrialMol& mol);
      double GenerateDihedrals(double* angles, double* angleEnergy,
			       double* angleWeights);
      void UseOldDih(double& energy, double& weight);

      DCData* data;
      uint atom, focus, prev, prevprev;
      uint angleKind, dihKind;
      double bondLength;
      double theta, phi;
      double bendEnergy, bendWeight;
   };

}

#endif

