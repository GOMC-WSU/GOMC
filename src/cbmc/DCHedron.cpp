#define _USE_MATH_DEFINES
#include <math.h>
#include "DCHedron.h"
#include "TrialMol.h"
#include "DCData.h"
#include "MolSetup.h"
#include "Forcefield.h"
#include "PRNG.h"
#include "NumLib.h"
#include <numeric>
#include <cassert>

namespace
{
   struct FindA1
   {
      FindA1(uint x) : x(x) {};
      bool operator()(const mol_setup::Bond& b) { return (b.a1 == x); }
      uint x;
   };

   struct FindAngle
   {
      FindAngle(uint x, uint y) : x(x), y(y) {}
      uint y, x;
      bool operator()(const mol_setup::Angle& a)
      {
         return (a.a0 == x && a.a2 == y) || (a.a0 == y && a.a2 == x);
      }
   };

}

namespace cbmc
{


   DCHedron::DCHedron(DCData* data, const mol_setup::MolKind& kind,
		      uint focus, uint prev)
      : data(data), focus(focus), prev(prev)
   {
      using namespace mol_setup;
      using namespace std;
      vector<Bond> onFocus = AtomBonds(kind, focus);
      for(uint i = 0; i < onFocus.size(); ++i) 
      {
	    if (onFocus[i].a1 == prev) 
	    {
	       anchorKind = onFocus[i].kind;
               eqAnchorBond = data->ff.bonds.Length(onFocus[i].kind);
               break;
            }
      }
      
      onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
		    onFocus.end());
      vector<Bond> onPrev = AtomBonds(kind, prev);
      onPrev.erase(remove_if(onPrev.begin(), onPrev.end(), FindA1(focus)),
		   onPrev.end());
      nBonds = onFocus.size();
      for (uint i = 0; i < nBonds; ++i)
      {
         bonded[i] = onFocus[i].a1;
         eqBondLength[i] = data->ff.bonds.Length(onFocus[i].kind);
	 bondKinds[i] = onFocus[i].kind;
      }
      vector<Angle> angles = AtomMidAngles(kind, focus);
      for (uint i = 0; i < nBonds; ++i)
      {
         typedef vector<Angle>::const_iterator Aiter;
         Aiter free = find_if(angles.begin(), angles.end(),
			      FindAngle(prev, bonded[i]));
         assert(free != angles.end());
         angleKinds[i][i] = free->kind;
         for (uint j = i + 1; j < nBonds; ++j) {
            Aiter pair = find_if(angles.begin(), angles.end(),
				 FindAngle(bonded[i], bonded[j]));
            angleKinds[i][j] = pair->kind;
            angleKinds[j][i] = pair->kind;
         }
      }

      phi[0] = 0.0;
      phiWeight[0] = 1.0;
   }


   double DCHedron::GetWeight()
   {
      double result = 1;
      for(uint i = 0; i < nBonds; ++i)
      {
         result *= thetaWeight[i];
         result *= phiWeight[i];
      }
      return result;
   }

   //Randomly generate nTrials angles and save energy and weight
   void DCHedron::GenerateAnglesNew(TrialMol& newMol, uint molIndex,
				    uint kind, uint nTrials, uint bType)
   {
      double* nonbonded_1_3 =  data->nonbonded_1_3;
      uint i;
      double distSq, thetaFix;
      bool angleFix = false;
      std::fill_n(nonbonded_1_3, nTrials, 0.0);

      if(data->ff.angles->AngleFixed(kind))
      {
	angleFix = true;
	thetaFix = data->ff.angles->Angle(kind);
      }

      for (i = 0; i < nTrials; ++i)
      {
	if(angleFix)
	  data->angles[i] = thetaFix;
	else
	  data->angles[i] = data->prng.rand(M_PI);
      }

#ifdef _OPENMP 
#pragma omp parallel for default(shared) private(i, distSq)
#endif      
      for (i = 0; i < nTrials; ++i)
      {
	 data->angleEnergy[i] = data->ff.angles->Calc(kind, data->angles[i]);
	 distSq = newMol.AngleDist(newAnchorBond, newBondLength[bType],
				   data->angles[i]);
	 nonbonded_1_3[i] =
	   data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);

	 if(isnan(nonbonded_1_3[i]))
	    nonbonded_1_3[i] = num::BIGNUM;

         data->angleWeights[i] = exp((data->angleEnergy[i] + nonbonded_1_3[i])
				     * -data->ff.beta);
      }
   }

   void DCHedron::GenerateAnglesOld(TrialMol& oldMol, uint molIndex,
				    uint kind, uint nTrials, uint bType)
   {
      double* nonbonded_1_3 =  data->nonbonded_1_3;
      uint i;
      double distSq, thetaFix;
      bool angleFix = false;
      std::fill_n(nonbonded_1_3, nTrials, 0.0);
    
      if(data->ff.angles->AngleFixed(kind))
      {
	angleFix = true;
	thetaFix = data->ff.angles->Angle(kind);
      }

      for (i = 0; i < nTrials; ++i)
      {
	 if(angleFix)
	   data->angles[i] = thetaFix;
	 else
	   data->angles[i] = data->prng.rand(M_PI);
      }

#ifdef _OPENMP 
#pragma omp parallel for default(shared) private(i, distSq)
#endif  
      for (i = 0; i < nTrials; ++i)
      {
	 data->angleEnergy[i] = data->ff.angles->Calc(kind, data->angles[i]);
	 distSq = oldMol.AngleDist(oldAnchorBond, oldBondLength[bType],
				   data->angles[i]);
	 nonbonded_1_3[i] =
	   data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);

	 if(isnan(nonbonded_1_3[i]))
	    nonbonded_1_3[i] = num::BIGNUM;

         data->angleWeights[i] = exp((data->angleEnergy[i] + nonbonded_1_3[i])
				      * -data->ff.beta);
      }
   }

   void DCHedron::FreeAnglesNew(TrialMol& newMol, uint molIndex, uint nTrials)
   {
      for (uint i = 0; i < nBonds; ++i)
      {
	 GenerateAnglesNew(newMol, molIndex, angleKinds[i][i], nTrials, i);
         double stepWeight = std::accumulate(data->angleWeights,
					     data->angleWeights + nTrials, 
					     0.0);
         uint winner = data->prng.PickWeighted(data->angleWeights,
					       nTrials, stepWeight);
         theta[i] = data->angles[winner];
         bendEnergy += data->angleEnergy[winner];
	 oneThree += data->nonbonded_1_3[winner];
         thetaWeight[i] = stepWeight;
      }
   }

   void DCHedron::FreeAnglesOld(TrialMol& oldMol, uint molIndex, uint nTrials)
   {  
      for (uint i = 0; i < nBonds; ++i)
      {
	 GenerateAnglesOld(oldMol, molIndex, angleKinds[i][i], nTrials, i);
         double stepWeight = std::accumulate(data->angleWeights,
					     data->angleWeights + nTrials, 
					     0.0);
         thetaWeight[i] = stepWeight;
      }
   }

   void DCHedron::PrepareNew(TrialMol& newMol, uint molIndex)
   {
      
      bendEnergy = 0.0;
      oneThree = 0.0;
      SetNewBond(newMol);
      FreeAnglesNew(newMol, molIndex, data->nAngleTrials);
      ConstrainedAngles(newMol, molIndex, data->nAngleTrials);
   }

   void DCHedron::PrepareOld(TrialMol& oldMol, uint molIndex)
   {
      oneThree = 0.0;
      bendEnergy = 0.0;
      SetOldBond(oldMol);
      FreeAnglesOld(oldMol, molIndex, data->nAngleTrials - 1);
   }

   void DCHedron::SetNewBond(TrialMol& newMol)
   {
      newBondEnergy = 0.0;
      newBondWeight = 1.0;
      double tempEn, tempW;
      for (uint i = 0; i < nBonds; ++i)
      {
	 if(data->ff.bonds.bondFixed(bondKinds[i]))
	 {
	   newBondLength[i] = eqBondLength[i];
	   tempEn = data->ff.bonds.Calc(bondKinds[i], newBondLength[i]);
	   tempW = exp(-1 * data->ff.beta * tempEn);
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
	     
	     newBondLength[i] = bond * eqBondLength[i];
	     tempEn = data->ff.bonds.Calc(bondKinds[i], newBondLength[i]);
	     tempW = exp(-1 * data->ff.beta * tempEn);  
 
	   }while(tempW < data->prng.rand());
	 }
	 newBondEnergy += tempEn;
	 newBondWeight *= tempW;
      }
      //find the anchor
      if(newMol.AtomExists(focus) && newMol.AtomExists(prev))
      {
	 //if both prev and focus are already built, we just use them 
	 newAnchorBond = sqrt(newMol.GetDistSq(focus, prev));
	 //tempEn = data->ff.bonds.Calc(anchorKind, newAnchorBond);
	 //tempW = exp(-1 * data->ff.beta * tempEn); 
	 //oldBondEnergy += tempEn;
	 //oldBondWeight *= tempW;
      }
      else
      {
	 if(data->ff.bonds.bondFixed(anchorKind))
	 {
	    newAnchorBond = eqAnchorBond;
	    tempEn = data->ff.bonds.Calc(anchorKind, newAnchorBond);
	    tempW = exp(-1 * data->ff.beta * tempEn);
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
	       
	       newAnchorBond = bond * eqAnchorBond;
	       tempEn = data->ff.bonds.Calc(anchorKind, newAnchorBond);
	       tempW = exp(-1 * data->ff.beta * tempEn);  
	       
	    }while(tempW < data->prng.rand());
	 }
	 newBondEnergy += tempEn;
	 newBondWeight *= tempW;
      }
   }
   
   void DCHedron::SetOldBond(TrialMol& oldMol)
   {
      double tempEn, tempW;
      oldBondEnergy = 0.0;
      oldBondWeight = 1.0;
      //set bond distance for old molecule
      for (uint i = 0; i < nBonds; ++i)
      {
	oldBondLength[i] = sqrt(oldMol.GetDistSq(focus, bonded[i]));
	tempEn = data->ff.bonds.Calc(bondKinds[i], oldBondLength[i]);
	tempW = exp(-1 * data->ff.beta * tempEn);
	oldBondEnergy += tempEn;
	oldBondWeight *= tempW;
      }
      //add old anchor
      if(oldMol.AtomExists(focus) && oldMol.AtomExists(prev))
      {
	oldAnchorBond = sqrt(oldMol.GetDistSq(focus, prev));
      }
      else
      {
	oldAnchorBond = sqrt(oldMol.GetDistSq(focus, prev));
	tempEn = data->ff.bonds.Calc(anchorKind, oldAnchorBond);
	tempW = exp(-1 * data->ff.beta * tempEn); 
	oldBondEnergy += tempEn;
	oldBondWeight *= tempW;
      }
   }


   void DCHedron::IncorporateOld(TrialMol& oldMol, uint molIndex)
   {
      bendEnergy = 0.0;
      oneThree = 0.0;
      const Forcefield& ff = data->ff;
      for (uint b = 0; b < nBonds; ++b)
      {
	 
         oldMol.OldThetaAndPhi(bonded[b], focus, theta[b], phi[b]);
         double thetaEnergy = data->ff.angles->Calc(angleKinds[b][b], theta[b]);
	 double distSq = oldMol.GetDistSq(prev, bonded[b]);
	 double nonbondedEn = data->calc.IntraEnergy_1_3(distSq, prev,
							 bonded[b], molIndex);

         thetaWeight[b] += exp(-1* data->ff.beta * (thetaEnergy + nonbondedEn));
         bendEnergy += thetaEnergy;
	 oneThree += nonbondedEn;
	 
	 if (b!=0)
	 {
	    double phiEnergy = 0.0;
	    nonbondedEn = 0.0;
	    phiWeight[b] = 0.0;
	    for (uint c = 0; c < b; ++c)
	    {
	       double cosTerm = cos(theta[b]) * cos(theta[c]);
	       double sinTerm = sin(theta[b]) * sin(theta[c]);
	       double bfcTheta = acos(sinTerm * cos(phi[b] - phi[c]) + 
				      cosTerm);
	
	       double distSq = oldMol.GetDistSq(bonded[c], bonded[b]);
	       nonbondedEn +=  data->calc.IntraEnergy_1_3(distSq, bonded[c],
						      bonded[b], molIndex);
	       
	       phiEnergy += ff.angles->Calc(angleKinds[b][c], bfcTheta);
	       
	    }
	    phiWeight[b] = exp(-ff.beta * (phiEnergy + nonbondedEn));
	    bendEnergy += phiEnergy;
	    oneThree += nonbondedEn;
	 }
      }
   }


      void DCHedron::ConstrainedAngles(TrialMol& newMol, uint molIndex,
				       uint nTrials)
   {
      double* angles = data->angles;
      double* energies = data->angleEnergy;
      double* weights = data->angleWeights;
      double* nonbonded_1_3 =  data->nonbonded_1_3;
      std::fill_n(nonbonded_1_3, nTrials, 0.0);
      phi[0] = 0.0;

      for (uint b = 1; b < nBonds; ++b)
      {
         //pick "twist" angles 
         for (uint i = 0; i < nTrials; ++i)
	 {
	    angles[i] = data->prng.rand(M_PI * 2); 
	    energies[i] = 0.0;
	    nonbonded_1_3[i] = 0.0;
         }

         //compare to angles determined in previous iterations
         for (uint c = 0; c < b; ++c)
	 {
            double cosTerm = cos(theta[b]) * cos(theta[c]);
            double sinTerm = sin(theta[b]) * sin(theta[c]);  
	    
            for (uint i = 0; i < nTrials; ++i)
            {
	       if(!data->ff.angles->AngleFixed(angleKinds[b][c]))
	       {
		 double bfcTheta = acos(sinTerm * cos(angles[i] - phi[c]) +
					cosTerm);
		 double distSq = newMol.AngleDist(newBondLength[b],
						  newBondLength[c], bfcTheta);
		 double tempEn = data->calc.IntraEnergy_1_3(distSq, bonded[b],
							    bonded[c],
							    molIndex);

		 if(isnan(tempEn))
		   tempEn = num::BIGNUM;

		 nonbonded_1_3[i] += tempEn;

		 energies[i] += data->ff.angles->Calc(angleKinds[b][c],
						      bfcTheta);
	       }
	       else
	       {
		  double fixedbfc = data->ff.angles->Angle(angleKinds[b][c]);
		  angles[i] = acos((cos(fixedbfc)-abs(cosTerm))/sinTerm)+phi[c];
		  double bfcTheta = acos(sinTerm * cos(angles[i] - phi[c]) +
					cosTerm);
		  double distSq = newMol.AngleDist(newBondLength[b],
						   newBondLength[c], bfcTheta);
		  double tempEn = data->calc.IntraEnergy_1_3(distSq, bonded[b],
							     bonded[c],
							     molIndex);
		  if(isnan(tempEn))
		    tempEn = num::BIGNUM;

		  nonbonded_1_3[i] += tempEn;
		  energies[i] += data->ff.angles->Calc(angleKinds[b][c],
						       bfcTheta);
	       }
            }
         }

         //calculate weights from combined energy
         double stepWeight = 0.0;
	 uint i;
#ifdef _OPENMP 
#pragma omp parallel for default(shared) private(i) reduction(+:stepWeight)
#endif 
         for (i = 0; i < nTrials; ++i)
	 {
	   weights[i] = exp(-1 * data->ff.beta * (energies[i] +
						  nonbonded_1_3[i]));
            stepWeight += weights[i];
         }

         uint winner = data->prng.PickWeighted(weights, nTrials, stepWeight);
         phi[b] = angles[winner];
         bendEnergy += energies[winner];
	 oneThree += nonbonded_1_3[winner];
         phiWeight[b] = stepWeight;
      }
   }

   //Calculate OldMol Bond Energy &
   //Calculate phi weight for nTrials using actual theta of OldMol
   void DCHedron::ConstrainedAnglesOld(uint nTrials, TrialMol& oldMol,
				       uint molIndex)
   {
      IncorporateOld(oldMol, molIndex);

      for (uint b = 1; b < nBonds; ++b) 
      {
	 double stepWeight = 0.0;
	 //pick "twist" angles 
	 for (uint i = 0; i < nTrials; ++i) 
	 {
	    double angles  = data->prng.rand(M_PI * 2);
	    double energies = 0.0;
	    double nonbondedEng = 0.0;
	    //compare to angles determined in previous iterations
	    for (uint c = 0; c < b; ++c) 
	    {	       
	       if(!data->ff.angles->AngleFixed(angleKinds[b][c]))
	       {
		 double cosTerm = cos(theta[b]) * cos(theta[c]);
		 double sinTerm = sin(theta[b]) * sin(theta[c]);
		 double bfcTheta = acos(sinTerm * cos(angles - phi[c]) 
					+ cosTerm);
		 double distSq = oldMol.AngleDist(oldBondLength[b],
						oldBondLength[c], bfcTheta);
		 nonbondedEng += data->calc.IntraEnergy_1_3(distSq, bonded[b],
							   bonded[c], molIndex);
		 if(isnan(nonbondedEng))
		   nonbondedEng = num::BIGNUM;

		 energies += data->ff.angles->Calc(angleKinds[b][c], bfcTheta);
	       }
	       else
	       {
		   double cosTerm = cos(theta[b]) * cos(theta[c]);
		   double sinTerm = sin(theta[b]) * sin(theta[c]);
		   double fixedbfc = data->ff.angles->Angle(angleKinds[b][c]);
		   angles = acos((cos(fixedbfc)-abs(cosTerm))/sinTerm)+phi[c];
		   double bfcTheta = acos(sinTerm * cos(angles - phi[c]) 
					  + cosTerm);
		   double distSq = oldMol.AngleDist(oldBondLength[b],
						    oldBondLength[c], bfcTheta);
		   nonbondedEng += data->calc.IntraEnergy_1_3(distSq, bonded[b],
							   bonded[c], molIndex);
		   if(isnan(nonbondedEng))
		     nonbondedEng = num::BIGNUM;

		   energies += data->ff.angles->Calc(angleKinds[b][c],
						     bfcTheta);
	       }
	    }
	    
	    //calculate weights from combined energy
	    double weights = exp(-1 * data->ff.beta *(energies + nonbondedEng));
	    stepWeight += weights;
	 }
	 phiWeight[b] += stepWeight;
      }
   }
}
