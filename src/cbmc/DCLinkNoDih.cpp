/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCLinkNoDih.h"
#include "TrialMol.h"
#include "../Forcefield.h"
#include "../XYZArray.h"
#include "../MoleculeKind.h"
#include "../MolSetup.h"


namespace cbmc
{
   DCLinkNoDih::DCLinkNoDih(DCData* data, const mol_setup::MolKind kind,
			    uint atom, uint focus)
     : data(data), atom(atom), focus(focus)
   {
      using namespace mol_setup;
      std::vector<Bond> bonds = AtomBonds(kind, atom);
      for (uint i = 0; i < bonds.size(); ++i)
      {
         if (bonds[i].a0 == focus || bonds[i].a1 == focus)
	 {
            bondLength = data->ff.bonds.Length(bonds[i].kind);
            break;
         }
      }
      std::vector<Angle> angles = AtomEndAngles(kind, atom);
      for (uint i = 0; i < angles.size(); ++i)
      {
         if (angles[i].a1 == focus)
	 {
            prev = angles[i].a2;
            angleKind = angles[i].kind;
			if (data->ff.angles->AngleEnergy(angleKind) > 999999999){
				angleFixe = true;
				thetaFix = data->ff.angles->Angle(angleKind);
			}
            break;
         }
      }
   }

   void DCLinkNoDih::PrepareNew()
   {
      double* angles = data->angles;
      double* angleEnergy = data->angleEnergy;
      double* angleWeights = data->angleWeights;
      PRNG& prng = data->prng;
      const Forcefield& ff = data->ff;
      uint count = data->nAngleTrials;
      bendWeight = 0;
      for (uint trial = 0; trial < count; trial++)
      {
		  if (angleFixe){
			  angles[trial] = thetaFix;
		  }
		  else{
			  angles[trial] = prng.rand(M_PI);
		  }
         angleEnergy[trial] = ff.angles->Calc(angleKind, angles[trial]);
         angleWeights[trial] = exp(angleEnergy[trial] * -ff.beta);
         bendWeight += angleWeights[trial];
      }
      uint winner = prng.PickWeighted(angleWeights, count, bendWeight);
      theta = angles[winner];
      bendEnergy = angleEnergy[winner];
   }

   void DCLinkNoDih::PrepareOld()
   {
      PRNG& prng = data->prng;
      const Forcefield& ff = data->ff;
      uint count = data->nAngleTrials - 1;
      bendWeight = 0;
      for (uint trial = 0; trial < count; trial++)
      {
		 double trialAngle;
		 if(angleFixe){
			 trialAngle = thetaFix;
		 }
		 else{
			 trialAngle = prng.rand(M_PI);
		 }
         double trialEn = ff.angles->Calc(angleKind, trialAngle);
         double trialWeight = exp(-ff.beta * trialEn);
         bendWeight += trialWeight;
      }
   }

   void DCLinkNoDih::IncorporateOld(TrialMol& oldMol)
   {
      double dummy;
      oldMol.OldThetaAndPhi(atom, focus, theta, dummy);
      const Forcefield& ff = data->ff;
      bendEnergy = ff.angles->Calc(angleKind, theta);
      bendWeight += exp(-ff.beta * bendEnergy);
   }

   void DCLinkNoDih::AlignBasis(TrialMol& mol)
   {
      mol.SetBasis(focus, prev);
   }

   void DCLinkNoDih::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      AlignBasis(oldMol);
      IncorporateOld(oldMol);
      double* nonbonded_1_4 = data->nonbonded_1_4;
      double* inter = data->inter;
	  double* real = data->real;
	  double *self = data->self;
	  double* corr = data->correction;
      uint nLJTrials = data->nLJTrialsNth;
      XYZArray& positions = data->positions;
      PRNG& prng = data->prng;
      positions.Set(0, oldMol.AtomPosition(atom));
      for (uint trial = 1, count = nLJTrials; trial < count; ++trial)
      {
         double phi = prng.rand(M_PI * 2);
         positions.Set(trial, oldMol.GetRectCoords(bondLength, theta, phi));
      }

      data->axes.WrapPBC(positions, oldMol.GetBox());
      std::fill_n(inter, nLJTrials, 0);
	  std::fill_n(self, nLJTrials, 0.0);
	  std::fill_n(real, nLJTrials, 0.0);
	  std::fill_n(corr, nLJTrials, 0.0);
      std::fill_n(nonbonded_1_4, nLJTrials, 0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
      data->calc.ParticleNonbonded_1_4(nonbonded_1_4, oldMol, positions, atom,
				   oldMol.GetBox(), nLJTrials);
	  data->calc.SwapSelf(self, molIndex, atom, oldMol.GetBox(), nLJTrials);
	  data->calc.SwapCorrection(corr, oldMol, positions, atom, 
					oldMol.GetBox(), nLJTrials);
      double stepWeight = 0;
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial)
      {
	stepWeight += exp(-data->ff.beta * (inter[trial] + nonbonded_1_4[trial]
				+ real[trial] + self[trial] + corr[trial]));
      }
      oldMol.MultWeight(stepWeight * bendWeight);
      oldMol.ConfirmOldAtom(atom);
      oldMol.AddEnergy(Energy(bendEnergy, nonbonded_1_4[0], inter[0],
				real[0], 0.0, self[0], corr[0]));
   }

   void DCLinkNoDih::BuildNew(TrialMol& newMol, uint molIndex)
   {
      AlignBasis(newMol);
      double* ljWeights = data->ljWeights;
      double* nonbonded_1_4 = data->nonbonded_1_4;
      double* inter = data->inter;
	  double *real = data->real;
	  double *self = data->self;
	  double* corr = data->correction;
      uint nLJTrials = data->nLJTrialsNth;
      XYZArray& positions = data->positions;
      PRNG& prng = data->prng;

      for (uint trial = 0, count = nLJTrials; trial < count; ++trial)
      {
         double phi = prng.rand(M_PI * 2);
         positions.Set(trial, newMol.GetRectCoords(bondLength, theta, phi));
      }

      data->axes.WrapPBC(positions, newMol.GetBox());
      std::fill_n(inter, nLJTrials, 0);
	  std::fill_n(self, nLJTrials, 0.0);
	  std::fill_n(real, nLJTrials, 0.0);
	  std::fill_n(corr, nLJTrials, 0.0);
	  std::fill_n(nonbonded_1_4, nLJTrials, 0);
	  std::fill_n(ljWeights, nLJTrials, 0.0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               newMol.GetBox(), nLJTrials);
      data->calc.ParticleNonbonded_1_4(nonbonded_1_4, newMol, positions, atom,
				   newMol.GetBox(), nLJTrials);
	  data->calc.SwapSelf(self, molIndex, atom, newMol.GetBox(), nLJTrials);
	  data->calc.SwapCorrection(corr, newMol, positions, atom, 
					newMol.GetBox(), nLJTrials);
      double stepWeight = 0.0;
      double beta = data->ff.beta;
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial)
      {
	 ljWeights[trial] = exp(-data->ff.beta * (inter[trial] + nonbonded_1_4[trial]
		 + real[trial] + self[trial] + corr[trial]));
         stepWeight += ljWeights[trial];
      }

      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      newMol.MultWeight(stepWeight * bendWeight);
      newMol.AddAtom(atom, positions[winner]);
      newMol.AddEnergy(Energy(bendEnergy, nonbonded_1_4[winner],
		  inter[winner], real[winner], 0.0, self[winner], corr[winner]));
   }

}

