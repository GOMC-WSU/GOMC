#include "DCOnSphere.h" 
#include "TrialMol.h" 
#include "DCData.h" 
#include "XYZArray.h" 
#include "PRNG.h" 
#include "Forcefield.h" 
#include "MolSetup.h" 
#include <omp.h>
 
namespace cbmc 
{ 
   DCOnSphere::DCOnSphere(DCData* data, const mol_setup::MolKind kind, 
			  uint atom, uint focus) : 
     data(data), atom(atom), 
     focus(focus), bondFix(true) 
   {  
      using namespace mol_setup; 
      std::vector<Bond> bonds = AtomBonds(kind, atom); 
      for(uint i = 0; i < bonds.size(); ++i) 
      { 
         if(bonds[i].a0 == focus || bonds[i].a1 == focus) 
	 { 
            eqBondLength = data->ff.bonds.Length(bonds[i].kind); 
	    bondKind = bonds[i].kind;
	    bondFix = data->ff.bonds.bondFixed(bondKind);
            break; 
         } 
      } 
   } 
 
   void DCOnSphere::SetOldBond(TrialMol& oldMol) 
   { 
     double BondDistSq = oldMol.OldDistSq(focus, atom); 
     oldBondLength = sqrt(BondDistSq);
     oldBondEnergy = data->ff.bonds.Calc(bondKind, oldBondLength);
     oldBondWeight = exp(-1 * data->ff.beta * oldBondEnergy);
   } 

   void DCOnSphere::SetNewBond(TrialMol& newMol)
   {
     if(bondFix)
     {
       newBondLength = eqBondLength;
       newBondEnergy = data->ff.bonds.Calc(bondKind, newBondLength);
       newBondWeight = exp(-1 * data->ff.beta * newBondEnergy);
     }
     else
     {
        double bond, bf;
        do
        {
	   do
	   {
	      bond = 0.2 * data->prng.rand() + 0.9;
	      bf = bond * bond * bond / 1.331;

	   }while(bf < data->prng.rand());

	   newBondLength = bond * eqBondLength;
	   newBondEnergy = data->ff.bonds.Calc(bondKind, newBondLength);
	   newBondWeight = exp(-1 * data->ff.beta * newBondEnergy);  
 
	}while(newBondWeight < data->prng.rand());
     }
   }
 
   void DCOnSphere::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* inter = data->inter; 
      double* self = data->self; 
      double* real = data->real; 
      double* correction = data->correction; 
      double stepWeight = 0; 
 
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(self, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(correction, nLJTrials, 0.0); 
      //considering bond energy for old molecule. There is no need to calculate 
      //for new molecule since we dont sample bond. 
      SetOldBond(oldMol); 
 
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, oldBondLength, 
					oldMol.AtomPosition(focus)); 
      positions.Set(0, oldMol.AtomPosition(atom)); 
      data->axes.WrapPBC(positions, oldMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               oldMol.GetBox(), nLJTrials);

#ifdef _OPENMP
#pragma omp parallel sections
#endif
{ 
#ifdef _OPENMP   
#pragma omp section
#endif
      data->calcEwald->SwapSelf(self, molIndex, atom, oldMol.GetBox(), 
			       nLJTrials);
#ifdef _OPENMP
#pragma omp section
#endif
      data->calcEwald->SwapCorrection(correction, oldMol, positions, atom,  
				     oldMol.GetBox(), nLJTrials); 
}
       
      const MoleculeKind& thisKind = oldMol.GetKind(); 
      double tempEn = 0.0; 
      for (uint i = 0; i < thisKind.NumAtoms(); i++) 
      { 
	 if (oldMol.AtomExists(i) && i != atom) 
	 { 
	    double distSq = oldMol.OldDistSq(i, atom); 
	    tempEn += data->calcEwald->CorrectionOldMol(oldMol, distSq, 
							     i, atom); 
	 } 
      } 
      correction[0] = tempEn; 
 
      for (uint trial = 0; trial < nLJTrials; trial++) 
      { 
         stepWeight += exp(-1 * data->ff.beta * 
			   (inter[trial] + real[trial] + self[trial] + 
			    correction[trial])); 
      } 
      oldMol.MultWeight(stepWeight * oldBondWeight); 
      oldMol.AddEnergy(Energy(oldBondEnergy, 0.0, inter[0], real[0], 0.0, 
			      self[0], correction[0])); 
      oldMol.ConfirmOldAtom(atom); 
   } 
 
   void DCOnSphere::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* inter = data->inter; 
      double* self = data->self; 
      double* real = data->real; 
      double* correction = data->correction; 
      double* ljWeights = data->ljWeights; 
      double stepWeight = 0; 
       
      std::fill_n(inter, nLJTrials, 0.0); 
      std::fill_n(self, nLJTrials, 0.0); 
      std::fill_n(real, nLJTrials, 0.0); 
      std::fill_n(correction, nLJTrials, 0.0); 
      std::fill_n(ljWeights, nLJTrials, 0.0); 

      SetNewBond(newMol);
 
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, newBondLength, 
					newMol.AtomPosition(focus)); 
      data->axes.WrapPBC(positions, newMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               newMol.GetBox(), nLJTrials);
#ifdef _OPENMP
#pragma omp parallel sections
#endif
{   
#ifdef _OPENMP   
#pragma omp section
#endif
      data->calcEwald->SwapSelf(self, molIndex, atom, newMol.GetBox(), 
			       nLJTrials); 
#ifdef _OPENMP
#pragma omp section
#endif
      data->calcEwald->SwapCorrection(correction, newMol, positions, atom,  
				     newMol.GetBox(), nLJTrials); 
}

      for (uint trial = 0; trial < nLJTrials; trial++) 
      { 
	 ljWeights[trial] = exp(-1 * data->ff.beta * 
			       (inter[trial] + real[trial] + self[trial] + 
				correction[trial])); 
         stepWeight += ljWeights[trial]; 
      } 
      uint winner = data->prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 
      newMol.MultWeight(stepWeight * newBondWeight); 
      newMol.AddEnergy(Energy(newBondEnergy, 0, inter[winner], real[winner],
			      0.0, self[winner], correction[winner])); 
      newMol.AddAtom(atom, positions[winner]); 
   }    
}      
