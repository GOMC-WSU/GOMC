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
	  double* self = data->self;
	  double *real = data->real;
	  double* corr = data->correction;
      std::fill_n(inter, nLJTrials, 0.0);
	  std::fill_n(self, nLJTrials, 0.0);
	  std::fill_n(real, nLJTrials, 0.0);
	  std::fill_n(corr, nLJTrials, 0.0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
	  if(DoEwald)
		data->calc.SwapSelf(self, molIndex, atom, oldMol.GetBox(), nLJTrials);
	  
	  double stepWeight = 0.0;
	  for(uint trial = 0; trial < nLJTrials; ++trial)
      {
	stepWeight += exp(-1 * data->ff.beta * (inter[trial] + real[trial] + self[trial]) );
      }
      oldMol.MultWeight(stepWeight);
      oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 0.0, self[0], 0.0));
      oldMol.ConfirmOldAtom(atom);
   }

   void DCSingle::BuildNew(TrialMol& newMol, uint molIndex)
   {
      PRNG& prng = data->prng;
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsFirst;
      double* inter = data->inter;
      double* ljWeights = data->ljWeights;
	  double* self = data->self;
	  double *real = data->real;
	  double *corr = data->correction;

      prng.FillWithRandom(positions, nLJTrials,
			  data->axes.GetAxis(newMol.GetBox()));
      std::fill_n(inter, nLJTrials, 0.0);
	  std::fill_n(self, nLJTrials, 0.0);
	  std::fill_n(real, nLJTrials, 0.0);
	  std::fill_n(corr, nLJTrials, 0.0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               newMol.GetBox(), nLJTrials);
	  if(DoEwald){
		data->calc.SwapSelf(self, molIndex, atom, newMol.GetBox(), nLJTrials);
	  }
      double stepWeight = 0.0;
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
	ljWeights[trial] = exp(-1 * data->ff.beta * (inter[trial] + real[trial] + self[trial]) );
	stepWeight += ljWeights[trial];
      }
      
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      double WinEnergy = inter[winner]+real[winner]+self[winner];
      newMol.MultWeight(stepWeight);
      newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner], 0.0, self[winner], 0.0));

      newMol.AddAtom(atom, positions[winner]);
   }
}
