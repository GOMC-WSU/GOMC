/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
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
         angles[trial] = prng.rand(M_PI);
         angleEnergy[trial] = ff.angles.Calc(angleKind, angles[trial]);
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
         double trialAngle = prng.rand(M_PI);
         double trialEn = ff.angles.Calc(angleKind, trialAngle);
         double trialWeight = exp(-ff.beta * trialEn);
         bendWeight += trialWeight;
      }
   }

   void DCLinkNoDih::IncorporateOld(TrialMol& oldMol)
   {
      double dummy;
      oldMol.OldThetaAndPhi(atom, focus, theta, dummy);
      const Forcefield& ff = data->ff;
      bendEnergy = ff.angles.Calc(angleKind, theta);
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
      double* inter = data->inter;
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
      data->calc.ParticleInter(inter, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
      double stepWeight = 0;
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial)
      {
         stepWeight += exp(-data->ff.beta * inter[trial]);
      }
      oldMol.MultWeight(stepWeight * bendWeight);
      oldMol.ConfirmOldAtom(atom);
      oldMol.AddEnergy(Energy(bendEnergy, 0, inter[0]));
   }

   void DCLinkNoDih::BuildNew(TrialMol& newMol, uint molIndex)
   {  
      AlignBasis(newMol);
      double* ljWeights = data->ljWeights;
      double* inter = data->inter;
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
      data->calc.ParticleInter(inter, positions, atom, molIndex,
                               newMol.GetBox(), nLJTrials);

      double stepWeight = 0;
      double beta = data->ff.beta;
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial)
      {
         ljWeights[trial] = exp(-data->ff.beta * inter[trial]);
         stepWeight += ljWeights[trial];
      }

      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      newMol.MultWeight(stepWeight * bendWeight);
      newMol.AddAtom(atom, positions[winner]);
      newMol.AddEnergy(Energy(bendEnergy, 0, inter[winner]));
   }

}

