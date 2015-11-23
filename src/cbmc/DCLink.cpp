#include "DCLink.h"
#include "TrialMol.h"
#include "../Forcefield.h"
#include "../XYZArray.h"
#include "../MoleculeKind.h"
#include "../MolSetup.h"


namespace cbmc
{

   DCLink::DCLink(DCData* data, const mol_setup::MolKind kind,
		  uint atom, uint focus)
      : data(data), atom(atom), focus(focus)
   {
      //will fail quietly if not a part of a valid linear molecule,
      //but we checked that, right?
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
	    angleKind = angles[i].kind;
            break;
	 }
      }
      std::vector<Dihedral> dihs = AtomEndDihs(kind, atom);
      for (uint i = 0; i < dihs.size(); ++i)
      {
         if (dihs[i].a1 == focus)
	 {
            dihKind = dihs[i].kind;
            prev = dihs[i].a2;
            prevprev = dihs[i].a3;
            break;
         }
      }
   }

   void DCLink::PrepareNew()
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
         angleEnergy[trial] = ff.angles->Calc(angleKind, angles[trial]);
         angleWeights[trial] = exp(angleEnergy[trial] * -ff.beta);
         bendWeight += angleWeights[trial];
      }
      uint winner = prng.PickWeighted(angleWeights, count, bendWeight);
      theta = angles[winner];
      bendEnergy = angleEnergy[winner];
   }

   void DCLink::PrepareOld()
   {
      PRNG& prng = data->prng;
      const Forcefield& ff = data->ff;
      uint count = data->nAngleTrials - 1;
      bendWeight = 0;
      for (uint trial = 0; trial < count; trial++)
      {
         double trialAngle = prng.rand(M_PI);
         double trialEn = ff.angles->Calc(angleKind, trialAngle);
         double trialWeight = exp(-ff.beta * trialEn);
         bendWeight += trialWeight;
      }
   }

   void DCLink::IncorporateOld(TrialMol& oldMol)
   {
      oldMol.OldThetaAndPhi(atom, focus, theta, phi);
      const Forcefield& ff = data->ff;
      bendEnergy = ff.angles->Calc(angleKind, theta);
      bendWeight += exp(-ff.beta * bendEnergy);
   }

   void DCLink::AlignBasis(TrialMol& mol)
   {
      mol.SetBasis(focus, prev, prevprev);
   }

   void DCLink::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      AlignBasis(oldMol);
      IncorporateOld(oldMol);
      double* angles = data->angles;
      double* angleEnergy = data->angleEnergy;
      double* angleWeights = data->angleWeights;
      double* ljWeights = data->ljWeights;
      double* torsion = data->bonded;
      double* nonbonded = data->nonbonded;
      double* inter = data->inter;
	  double *real = data->real;
	  double *self = data->self;
	  double* corr = data->correction;
      uint nLJTrials = data->nLJTrialsNth;
      XYZArray& positions = data->positions;
      PRNG& prng = data->prng;

      UseOldDih(torsion[0], ljWeights[0]);
      positions.Set(0, oldMol.AtomPosition(atom));
      for (uint trial = 1; trial < nLJTrials; ++trial)
      {
         ljWeights[trial] = GenerateDihedrals(angles, angleEnergy,
					      angleWeights);
         uint winner = prng.PickWeighted(angleWeights, data->nDihTrials,
					ljWeights[trial]);
         torsion[trial] = angleEnergy[winner];
         positions.Set(trial, oldMol.GetRectCoords(bondLength, theta,
						   angles[winner]));
      }
      data->axes.WrapPBC(positions, oldMol.GetBox());
      std::fill_n(inter, nLJTrials, 0);
	  std::fill_n(self, nLJTrials, 0);
	  std::fill_n(real, nLJTrials, 0);
      std::fill_n(nonbonded, nLJTrials, 0);
	  std::fill_n(corr, nLJTrials, 0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               oldMol.GetBox(), nLJTrials);
      data->calc.ParticleNonbonded(nonbonded, oldMol, positions, atom,
				   oldMol.GetBox(), nLJTrials);
	  if(DoEwald){
		data->calc.SwapSelf(self, molIndex, atom, oldMol.GetBox(), nLJTrials);
		data->calc.SwapCorrection(corr, oldMol, positions, atom, oldMol.GetBox(), nLJTrials);
	  }

      double dihLJWeight = 0;
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         ljWeights[trial] *= exp(-data->ff.beta *
				 (inter[trial] + nonbonded[trial] + real[trial] + self[trial] + corr[trial]));
         dihLJWeight += ljWeights[trial];
      }
      oldMol.MultWeight(dihLJWeight * bendWeight);
      oldMol.ConfirmOldAtom(atom);
      oldMol.AddEnergy(Energy(torsion[0] + bendEnergy, nonbonded[0],
			      inter[0], real[0], 0, self[0], corr[0]));
   }

   void DCLink::BuildNew(TrialMol& newMol, uint molIndex)
   {
      AlignBasis(newMol);
      double* angles = data->angles;
      double* angleEnergy = data->angleEnergy;
      double* angleWeights = data->angleWeights;
      double* ljWeights = data->ljWeights;
      double* torsion = data->bonded;
      double* nonbonded = data->nonbonded;
      double* inter = data->inter;
	  double *real = data->real;
	  double* self = data->self;
	  double* corr = data->correction;
      uint nLJTrials = data->nLJTrialsNth;
      XYZArray& positions = data->positions;
      PRNG& prng = data->prng;

      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         ljWeights[trial] = GenerateDihedrals(angles, angleEnergy,
					      angleWeights);
         uint winner = prng.PickWeighted(angleWeights, data->nDihTrials,
					 ljWeights[trial]);
         torsion[trial] = angleEnergy[winner];
         positions.Set(trial, newMol.GetRectCoords(bondLength, theta, 
						   angles[winner]));
      }
      data->axes.WrapPBC(positions, newMol.GetBox());
      std::fill_n(inter, nLJTrials, 0);
      std::fill_n(nonbonded, nLJTrials, 0);
	  std::fill_n(self, nLJTrials, 0);
	  std::fill_n(real, nLJTrials, 0);
	  std::fill_n(corr, nLJTrials, 0);
      data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                               newMol.GetBox(), nLJTrials);
      data->calc.ParticleNonbonded(nonbonded, newMol, positions, atom,
				   newMol.GetBox(), nLJTrials);
	  if(DoEwald){
		data->calc.SwapSelf(self, molIndex, atom, newMol.GetBox(), nLJTrials);
		data->calc.SwapCorrection(corr, newMol, positions, atom, newMol.GetBox(), nLJTrials);
	  }

      double dihLJWeight = 0;
      double beta = data->ff.beta;
      for (uint trial = 0; trial < nLJTrials; ++trial)
      {
         ljWeights[trial] *= exp(-data->ff.beta * (inter[trial] + nonbonded[trial] + real[trial] + self[trial] + corr[trial]));
         dihLJWeight += ljWeights[trial];
      }

      uint winner = prng.PickWeighted(ljWeights, nLJTrials, dihLJWeight);
      double WinEnergy = inter[winner]+real[winner]+self[winner]+corr[winner];
      
      newMol.MultWeight(dihLJWeight * bendWeight);
      newMol.AddAtom(atom, positions[winner]);
      newMol.AddEnergy(Energy(torsion[winner] + bendEnergy, nonbonded[winner],
			      inter[winner], real[winner], 0, self[winner], corr[winner]));
   }

   double DCLink::GenerateDihedrals(double* angles, double* angleEnergy,
				    double* angleWeights)
   {
      double stepWeight = 0.0;
      PRNG& prng = data->prng;
      const Forcefield& ff = data->ff;
      for (uint trial = 0, count = data->nDihTrials; trial < count; ++trial)
      {
         angles[trial] = prng.rand(2 * M_PI);
         angleEnergy[trial] = ff.dihedrals.Calc(dihKind, angles[trial]);
         angleWeights[trial] = exp(-ff.beta * angleEnergy[trial]);
         stepWeight += angleWeights[trial];
      }
      return stepWeight;
   }

   void DCLink::UseOldDih(double& energy, double& weight)
   {
      PRNG& prng = data->prng;
      const Forcefield& ff = data->ff;
      double beta = data->ff.beta;

      energy = ff.dihedrals.Calc(dihKind, phi);
      weight = exp(-beta * energy);
      for (uint trial = data->nDihTrials - 1; trial-- > 0;)
      {
         double trialPhi = prng.rand(2 * M_PI);
         double trialEnergy = ff.dihedrals.Calc(dihKind, trialPhi);
         weight += exp(-beta * trialEnergy);
      }
   }

}
