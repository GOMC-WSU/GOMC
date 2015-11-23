#include "DCOnSphere.h"
#include "TrialMol.h"
#include "DCData.h"
#include "../XYZArray.h"
#include "../PRNG.h"
#include "../Forcefield.h"
#include "../MolSetup.h"

namespace cbmc
{
   DCOnSphere::DCOnSphere(DCData* data, const mol_setup::MolKind kind,
			  uint atom, uint focus) :
     data(data), atom(atom),
     focus(focus)
   { 
      using namespace mol_setup;
      std::vector<Bond> bonds = AtomBonds(kind, atom);
      for(uint i = 0; i < bonds.size(); ++i)
      {
         if(bonds[i].a0 == focus || bonds[i].a1 == focus)
		{
            bondLength = data->ff.bonds.Length(bonds[i].kind);
            break;
         }
      }
   }

  void DCOnSphere::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsNth;
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
					oldMol.AtomPosition(focus));
      positions.Set(0, oldMol.AtomPosition(atom));

      double* inter = data->inter;
      double* self = data->self;
      double *real = data->real;
      double* corr = data->correction;
      double stepWeight = 0.0;
      data->axes.WrapPBC(positions, oldMol.GetBox());
      std::fill_n(inter, nLJTrials, 0.0);
	  std::fill_n(self, nLJTrials, 0.0);
	  std::fill_n(real, nLJTrials, 0.0);
	  std::fill_n(corr, nLJTrials, 0.0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
	  if(DoEwald){
		data->calc.SwapSelf(self, molIndex, atom, oldMol.GetBox(), nLJTrials);
		data->calc.SwapCorrection(corr, oldMol, positions, atom, oldMol.GetBox(), nLJTrials);
	  }
      
      for (uint trial = 0; trial < nLJTrials; trial++)
      {
	stepWeight += exp(-1 * data->ff.beta * (inter[trial] + real[trial] + self[trial] + corr[trial]) );
      }
      oldMol.MultWeight(stepWeight);
      oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 0.0, self[0], corr[0]));
      oldMol.ConfirmOldAtom(atom);
   }

  void DCOnSphere::BuildNew(TrialMol& newMol, uint molIndex)
   {
      XYZArray& positions = data->positions;
      uint nLJTrials = data->nLJTrialsNth;
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
					newMol.AtomPosition(focus));

      double* inter = data->inter;
      double* self = data->self;
      double *real = data->real;
      double* corr = data->correction;
      double* ljWeights = data->ljWeights;
      double stepWeight = 0.0;
      data->axes.WrapPBC(positions, newMol.GetBox());
      
      std::fill_n(inter, nLJTrials, 0.0);
      std::fill_n(self, nLJTrials, 0.0);
      std::fill_n(real, nLJTrials, 0.0);
      std::fill_n(corr, nLJTrials, 0.0);
      std::fill_n(ljWeights, nLJTrials, 0.0);

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, newMol.GetBox(), nLJTrials);

      if(DoEwald){
	data->calc.SwapSelf(self, molIndex, atom, newMol.GetBox(), nLJTrials);
	data->calc.SwapCorrection(corr, newMol, positions, atom, newMol.GetBox(), nLJTrials);
      }
      for (uint trial = 0; trial < nLJTrials; trial++)
      {
         ljWeights[trial] = exp(-1 * data->ff.beta * (inter[trial] + real[trial] + self[trial] + corr[trial]) );
         stepWeight += ljWeights[trial];
      }
            
      uint winner = data->prng.PickWeighted(ljWeights, nLJTrials, stepWeight);

      double WinEnergy = inter[winner]+real[winner]+self[winner]+corr[winner];

      newMol.MultWeight(stepWeight);
      newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner], 0.0, self[winner], corr[winner]));
      newMol.AddAtom(atom, positions[winner]);
   }   
}
