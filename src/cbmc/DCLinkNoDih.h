/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version) 
Copyright (C) 2015  GOMC Group 
A copy of the GNU General Public License can be found in the COPYRIGHT.txt 
along with this program, also can be found at <http://www.gnu.org/licenses/>. 
*******************************************************************************/
#ifndef DCLINKNODIH_H 
#define DCLINKNODIH_H 
 
#include "BasicTypes.h" 
#include "DCComponent.h" 
#include "DCData.h" 
 
namespace cbmc 
{ 
   class DCLinkNoDih : public DCComponent 
   { 
    public: 
      DCLinkNoDih(DCData* data, const mol_setup::MolKind kind, 
		  uint atom, uint focus); 
 
      virtual void PrepareNew(TrialMol& newMol, uint molIndex); 
      virtual void PrepareOld(TrialMol& oldMol, uint molIndex); 
      virtual void BuildOld(TrialMol& oldMol, uint molIndex); 
      virtual void BuildNew(TrialMol& newMol, uint molIndex); 
 
      virtual DCLinkNoDih* Clone() { return new DCLinkNoDih(*this); }; 
 
    private: 
      void IncorporateOld(TrialMol& oldMol, uint molIndex); 
      void IncorporateNew(TrialMol& newMol, uint molIndex); 
      void AlignBasis(TrialMol& mol); 
      void SetOldBond(TrialMol& oldMol); 
      void SetNewBond(TrialMol& newMol);
 
      DCData* data; 

      uint atom, focus, prev; 
      uint bondKind, angleKind; 

      //need to habe for 1-3 interaction calculation
      double bond[2]; 
      double oldBond[2];
 
      double eqBondLength, oldBondLength, newBondLength; 
      double theta, thetaFix; 
      double bendEnergy, bendWeight, oneThree; 
      double oldBondEnergy, oldBondWeight;
      double newBondEnergy, newBondWeight;
      bool bondFix, angleFix;
   }; 
} 
 
#endif 
