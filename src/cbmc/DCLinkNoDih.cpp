#include "DCLinkNoDih.h" 
#include "TrialMol.h" 
#include "Forcefield.h" 
#include "XYZArray.h" 
#include "MoleculeKind.h" 
#include "MolSetup.h" 
#include "NumLib.h"
 
 
namespace cbmc 
{ 
   DCLinkNoDih::DCLinkNoDih(DCData* data, const mol_setup::MolKind kind, 
			    uint atom, uint focus) 
     : data(data), atom(atom), focus(focus), angleFix(false) 
   { 
      using namespace mol_setup; 
      std::vector<Bond> bonds = AtomBonds(kind, atom); 
      for (uint i = 0; i < bonds.size(); ++i) 
      { 
         if (bonds[i].a0 == focus || bonds[i].a1 == focus) 
	 { 
            bondLength = data->ff.bonds.Length(bonds[i].kind); 
	    bond[1] = bondLength; 
	    bondKind = bonds[i].kind; 
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
	    if (data->ff.angles->AngleFixed(angleKind)) 
	    { 
	       angleFix = true; 
	       thetaFix = data->ff.angles->Angle(angleKind); 
	    } 
            break; 
         } 
      } 
 
      std::vector<Bond> bond2 = AtomBonds(kind, prev); 
      for (uint i = 0; i < bond2.size(); ++i) 
      { 
         if (bond2[i].a0 == focus || bond2[i].a1 == focus) 
	 { 
            bond[0] = data->ff.bonds.Length(bond2[i].kind); 
            break; 
         } 
      } 
       
   } 
 
   void DCLinkNoDih::PrepareNew(TrialMol& newMol, uint molIndex) 
   { 
      double* angles = data->angles; 
      double* angleEnergy = data->angleEnergy; 
      double* angleWeights = data->angleWeights; 
      double* nonbonded_1_3 =  data->nonbonded_1_3; 
      PRNG& prng = data->prng; 
      const Forcefield& ff = data->ff; 
      uint count = data->nAngleTrials; 
      std::fill_n(nonbonded_1_3, count, 0.0); 
      bendWeight = 0.0; 
 
 
      for (uint trial = 0; trial < count; trial++) 
      { 
	 if (angleFix) 
	 { 
	    angles[trial] = thetaFix; 
	    angleEnergy[trial] = 0.0; 
	 } 
	 else 
	 { 
	    angles[trial] = prng.rand(M_PI); 
	    angleEnergy[trial] = ff.angles->Calc(angleKind, angles[trial]);  
	 } 

	 double distSq = newMol.AngleDist(bond[0], bond[1], angles[trial]); 
	 nonbonded_1_3[trial] = data->calc.IntraEnergy_1_3(distSq, prev, atom, 
							   molIndex);  
	 if(isnan(nonbonded_1_3[trial])) 
	   nonbonded_1_3[trial] = num::BIGNUM; 
	  
         angleWeights[trial] = exp((angleEnergy[trial] + nonbonded_1_3[trial]) 
				   * -ff.beta); 
         bendWeight += angleWeights[trial]; 
      } 
      uint winner = prng.PickWeighted(angleWeights, count, bendWeight); 
      theta = angles[winner]; 
      bendEnergy = angleEnergy[winner]; 
      oneThree = nonbonded_1_3[winner]; 
   } 
 
   void DCLinkNoDih::PrepareOld(TrialMol& oldMol, uint molIndex) 
   { 
      PRNG& prng = data->prng; 
      const Forcefield& ff = data->ff; 
      uint count = data->nAngleTrials - 1; 
      bendWeight = 0; 
       
      //set bond distance for old molecule 
      double BondDistSq1 = oldMol.OldDistSq(focus, atom); 
      double BondDistSq2 = oldMol.OldDistSq(prev, focus); 
      SetOldMolBond(1, BondDistSq1); 
      SetOldMolBond(0, BondDistSq2); 
 
      for (uint trial = 0; trial < count; trial++) 
      {
	 double trialAngle; 
	 double trialEn; 
	 if(angleFix) 
	 { 
	    trialAngle = thetaFix; 
	    trialEn = 0.0; 
	 } 
	 else 
	 { 
	    trialAngle = prng.rand(M_PI); 
	    trialEn = ff.angles->Calc(angleKind, trialAngle); 
	 } 
	 double distSq = oldMol.AngleDist(oldBond[0], oldBond[1], trialAngle); 
 
	 double tempEn = data->calc.IntraEnergy_1_3(distSq, prev, atom, molIndex); 
	 if(isnan(tempEn)) 
	   tempEn = num::BIGNUM;
         trialEn += tempEn; 
         double trialWeight = exp(-ff.beta * trialEn); 
         bendWeight += trialWeight; 
      } 
   } 
 
  void DCLinkNoDih::IncorporateOld(TrialMol& oldMol, uint molIndex) 
   { 
      double dummy; 
      oldMol.OldThetaAndPhi(atom, focus, theta, dummy); 
      const Forcefield& ff = data->ff; 

      bendEnergy = ff.angles->Calc(angleKind, theta); 
      double distSq = oldMol.OldDistSq(prev, atom); 
      oneThree = data->calc.IntraEnergy_1_3(distSq, prev, atom, molIndex); 
      bendWeight += exp(-ff.beta * (bendEnergy + oneThree)); 
      //considering bond energy for old molecule. There is no need to calculate 
      //for new molecule since we dont sample bond. 
      double BondDistSq = oldMol.OldDistSq(focus, atom); 
      oldBondEnergy = ff.bonds.Calc(bondKind, sqrt(BondDistSq)); 
   } 
 
   void DCLinkNoDih::AlignBasis(TrialMol& mol) 
   { 
      mol.SetBasis(focus, prev); 
   } 
 
   void DCLinkNoDih::SetOldMolBond(const uint i, const double distSq) 
   { 
     oldBond[i] = sqrt(distSq); 
   } 
 
   void DCLinkNoDih::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      AlignBasis(oldMol); 
      IncorporateOld(oldMol, molIndex); 
      double* inter = data->inter; 
      double* real = data->real;  
      uint nLJTrials = data->nLJTrialsNth; 
      XYZArray& positions = data->positions; 
      PRNG& prng = data->prng; 
 
      std::fill_n(inter, nLJTrials, 0.0);  
      std::fill_n(real, nLJTrials, 0.0);  
 
      positions.Set(0, oldMol.AtomPosition(atom)); 
      for (uint trial = 1, count = nLJTrials; trial < count; ++trial) 
      { 
         double phi = prng.rand(M_PI * 2); 
         positions.Set(trial, oldMol.GetRectCoords(bondLength, theta, phi)); 
      } 
 
      data->axes.WrapPBC(positions, oldMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               oldMol.GetBox(), nLJTrials); 
 
      double stepWeight = 0; 
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial) 
      { 
	stepWeight += exp(-data->ff.beta * (inter[trial] + real[trial])); 
      } 
      oldMol.MultWeight(stepWeight * bendWeight); 
      oldMol.ConfirmOldAtom(atom); 
      oldMol.AddEnergy(Energy(bendEnergy + oldBondEnergy, oneThree, inter[0], 
			      real[0], 0.0, 0.0, 0.0)); 
   } 
 
   void DCLinkNoDih::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      AlignBasis(newMol); 
      double* ljWeights = data->ljWeights; 
      double* inter = data->inter; 
      double* real = data->real; 
      uint nLJTrials = data->nLJTrialsNth; 
      XYZArray& positions = data->positions; 
      PRNG& prng = data->prng; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial) 
      { 
         double phi = prng.rand(M_PI * 2); 
         positions.Set(trial, newMol.GetRectCoords(bondLength, theta, phi)); 
      } 
 
      data->axes.WrapPBC(positions, newMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               newMol.GetBox(), nLJTrials);
 
      double stepWeight = 0; 
      double beta = data->ff.beta; 
      for (uint trial = 0, count = nLJTrials; trial < count; ++trial) 
      { 
	 ljWeights[trial] = exp(-data->ff.beta * 
				(inter[trial] + real[trial])); 
         stepWeight += ljWeights[trial]; 
      } 
 
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 
      newMol.MultWeight(stepWeight * bendWeight); 
      newMol.AddAtom(atom, positions[winner]); 
      newMol.AddEnergy(Energy(bendEnergy, oneThree, inter[winner], 
			      real[winner], 0.0, 0.0, 0.0)); 
   } 
 
}           
