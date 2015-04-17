/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCDiatomic.h"

#include "../MolSetup.h"
#include "DCData.h"
#include "../PRNG.h"
#include "TrialMol.h"
#include "../Forcefield.h"
#include "../CalculateEnergy.h"
#include "../EnergyTypes.h"
#include <cmath>

namespace cbmc{


   DCDiatomic::DCDiatomic(DCData* data, const mol_setup::MolKind kind,
			  uint first, uint second)
      : data(data), first(first), second(second)
   {
      using namespace mol_setup;
      std::vector<Bond> bonds = AtomBonds(kind, first);
      for(uint i = 0; i < bonds.size(); ++i) {
         if(bonds[i].a0 == second || bonds[i].a1 == second) {
            bondLength = data->ff.bonds.Length(bonds[i].kind);
            return;
         }
      }
   }

   void DCDiatomic::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      PRNG& prng = data->prng;
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsNth;
      double beta = data->ff.beta;

      prng.FillWithRandom(data->positions, nLJTrials,
			  data->axes.GetAxis(oldMol.GetBox()));
      positions.Set(0, oldMol.AtomPosition(first));

      double* inter = data->inter;
      double stepWeight = 0;
      data->axes.WrapPBC(positions, oldMol.GetBox());
      data->calc.ParticleInter(inter,positions,first,   molIndex,
			       oldMol.GetBox(), nLJTrials);
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         stepWeight += exp(-1 * beta * inter[trial]);
      }
      oldMol.MultWeight(stepWeight);
      oldMol.AddEnergy(Energy(0, 0, inter[0]));
      oldMol.ConfirmOldAtom(first);

      prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
				  oldMol.AtomPosition(first));
      positions.Set(0, oldMol.AtomPosition(second));

      stepWeight = 0;
      data->axes.WrapPBC(positions, oldMol.GetBox());
      data->calc.ParticleInter(inter, positions,second,  molIndex,
			       oldMol.GetBox(), nLJTrials);
      for (uint trial = 0; trial < nLJTrials; trial++)
      {
         stepWeight += exp(-1 * beta * inter[trial]);
      }
      oldMol.MultWeight(stepWeight);
      oldMol.AddEnergy(Energy(0, 0, inter[0]));
      oldMol.ConfirmOldAtom(second);
   }

   void DCDiatomic::BuildNew(TrialMol& newMol, uint molIndex)
   {   
      PRNG& prng = data->prng;
      XYZArray& positions = data->positions;
      double beta = data->ff.beta;
      uint nLJTrials = data->nLJTrialsNth;
      double* inter = data->inter;
      double* ljWeights = data->ljWeights;

      prng.FillWithRandom(positions, nLJTrials,
			  data->axes.GetAxis(newMol.GetBox()));
      data->axes.WrapPBC(positions, newMol.GetBox());
      data->calc.ParticleInter(inter, positions,first,  molIndex,
			       newMol.GetBox(), nLJTrials);

      double stepWeight = 0;
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         ljWeights[trial] = exp(-1 * beta * inter[trial]);
         stepWeight += ljWeights[trial];
      }
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      newMol.MultWeight(stepWeight);
      newMol.AddEnergy(Energy(0, 0, inter[winner]));
      newMol.AddAtom(first, positions[winner]);

      prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
				  newMol.AtomPosition(first));
      data->axes.WrapPBC(positions, newMol.GetBox());
      data->calc.ParticleInter(inter, positions,second,  molIndex,
			       newMol.GetBox(), nLJTrials);

      stepWeight = 0;
      for (uint trial = 0; trial < nLJTrials; trial++)
      {
         ljWeights[trial] = exp(-1 * beta * inter[trial]);
         stepWeight += ljWeights[trial];
      }
      winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      newMol.MultWeight(stepWeight);
      newMol.AddEnergy(Energy(0, 0, inter[winner]));
      newMol.AddAtom(second, positions[winner]);
   }

}

