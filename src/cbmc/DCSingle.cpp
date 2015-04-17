/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCSingle.h"
#include "TrialMol.h"
#include "DCData.h"
#include "../PRNG.h"
#include "../CalculateEnergy.h"
#include "../XYZArray.h"
#include "../Forcefield.h"

namespace cbmc
{

   void DCSingle::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      PRNG& prng = data->prng;
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsFirst;

      prng.FillWithRandom(data->positions, nLJTrials,
			  data->axes.GetAxis(oldMol.GetBox()));
      positions.Set(0, oldMol.AtomPosition(atom));

      double* inter = data->inter;
      double stepWeight = 0;
      std::fill_n(inter, nLJTrials, 0);
      data->calc.ParticleInter(inter, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         stepWeight += exp(-1 * data->ff.beta * inter[trial]);
      }
      oldMol.MultWeight(stepWeight);
      oldMol.AddEnergy(Energy(0, 0, inter[0]));
      oldMol.ConfirmOldAtom(atom);
   }

   void DCSingle::BuildNew(TrialMol& newMol, uint molIndex)
   { 
      PRNG& prng = data->prng;
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsFirst;
      double* inter = data->inter;
      double* ljWeights = data->ljWeights;

      prng.FillWithRandom(positions, nLJTrials,
			  data->axes.GetAxis(newMol.GetBox()));
      std::fill_n(inter, nLJTrials, 0);
      data->calc.ParticleInter(inter, positions, atom, molIndex,
                               newMol.GetBox(), nLJTrials);

      double stepWeight = 0;
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {   //printf("trail energy=%f\n", inter[trial]);
         ljWeights[trial] = exp(-1 * data->ff.beta * inter[trial]);
         stepWeight += ljWeights[trial];
      }
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
	  // printf("CPU trial winner=%f\n", inter[winner]);
      newMol.MultWeight(stepWeight);
      newMol.AddEnergy(Energy(0, 0, inter[winner]));
      newMol.AddAtom(atom, positions[winner]);
   }
}

