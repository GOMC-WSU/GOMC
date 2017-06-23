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
     focus(focus) 
   {  
      using namespace mol_setup; 
      std::vector<Bond> bonds = AtomBonds(kind, atom); 
      for(uint i = 0; i < bonds.size(); ++i) 
      { 
         if(bonds[i].a0 == focus || bonds[i].a1 == focus) 
	 { 
            bondLength = data->ff.bonds.Length(bonds[i].kind); 
	    bondKind = bonds[i].kind; 
            break; 
         } 
      } 
   } 
 
   void DCOnSphere::SetOldBondEnergy(TrialMol& oldMol) 
   { 
     double BondDistSq = oldMol.OldDistSq(focus, atom); 
     oldBondEnergy = data->ff.bonds.Calc(bondKind, sqrt(BondDistSq)); 
   } 
 
   void DCOnSphere::BuildOld(TrialMol& oldMol, uint molIndex) 
   { 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* inter = data->inter; 
      double* real = data->real; 
      double stepWeight = 0; 
 
      std::fill_n(inter, nLJTrials, 0.0);
      std::fill_n(real, nLJTrials, 0.0); 
      //considering bond energy for old molecule. There is no need to calculate 
      //for new molecule since we dont sample bond. 
      SetOldBondEnergy(oldMol); 
 
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength, 
					oldMol.AtomPosition(focus)); 
      positions.Set(0, oldMol.AtomPosition(atom)); 
      data->axes.WrapPBC(positions, oldMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               oldMol.GetBox(), nLJTrials);
 
 
      for (uint trial = 0; trial < nLJTrials; trial++) 
      { 
         stepWeight += exp(-1 * data->ff.beta * 
			   (inter[trial] + real[trial])); 
      } 
      oldMol.MultWeight(stepWeight); 
      oldMol.AddEnergy(Energy(oldBondEnergy, 0.0, inter[0], real[0], 0.0, 
			      0.0, 0.0)); 
      oldMol.ConfirmOldAtom(atom); 
   } 
 
   void DCOnSphere::BuildNew(TrialMol& newMol, uint molIndex) 
   { 
      XYZArray& positions = data->positions; 
      uint nLJTrials = data->nLJTrialsNth; 
      double* inter = data->inter; 
      double* real = data->real;  
      double* ljWeights = data->ljWeights; 
      double stepWeight = 0; 
       
      std::fill_n(inter, nLJTrials, 0.0);  
      std::fill_n(real, nLJTrials, 0.0);  
      std::fill_n(ljWeights, nLJTrials, 0.0); 
 
      data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength, 
					newMol.AtomPosition(focus)); 
      data->axes.WrapPBC(positions, newMol.GetBox()); 

      data->calc.ParticleInter(inter, real, positions, atom, molIndex, 
                               newMol.GetBox(), nLJTrials);


      for (uint trial = 0; trial < nLJTrials; trial++) 
      { 
	 ljWeights[trial] = exp(-1 * data->ff.beta * 
			       (inter[trial] + real[trial])); 
         stepWeight += ljWeights[trial]; 
      } 
      uint winner = data->prng.PickWeighted(ljWeights, nLJTrials, stepWeight); 
      newMol.MultWeight(stepWeight); 
      newMol.AddEnergy(Energy(0, 0, inter[winner], real[winner], 0.0,
			      0.0, 0.0)); 
      newMol.AddAtom(atom, positions[winner]); 
   }    
}      
