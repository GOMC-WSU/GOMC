/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#define _USE_MATH_DEFINES
#include <math.h>
#include "DCHedron.h"
#include "TrialMol.h"
#include "DCData.h"
#include "../MolSetup.h"
#include "../Forcefield.h"
#include "../PRNG.h"
#include <numeric>
#include <cassert>

namespace
{
   //Wish I could use lambdas..
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
      onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
		    onFocus.end());
      vector<Bond> onPrev = AtomBonds(kind, prev);
      onPrev.erase(remove_if(onPrev.begin(), onPrev.end(), FindA1(focus)),
		   onPrev.end());
      nBonds = onFocus.size();
      for (uint i = 0; i < nBonds; ++i)
      {
         bonded[i] = onFocus[i].a1;
         bondLength[i] = data->ff.bonds.Length(onFocus[i].kind);
      }
      vector<Angle> angles = AtomMidAngles(kind, focus);
      for (uint i = 0; i < nBonds; ++i)
      {
         typedef vector<Angle>::const_iterator Aiter;
         Aiter free = find_if(angles.begin(), angles.end(), FindAngle(prev, bonded[i]));
         assert(free != angles.end());
         angleKinds[i][i] = free->kind;
         for (uint j = i + 1; j < nBonds; ++j) {
            Aiter pair = find_if(angles.begin(), angles.end(), FindAngle(bonded[i], bonded[j]));
            angleKinds[i][j] = pair->kind;
            angleKinds[j][i] = pair->kind;
         }
      }

      phi[0] = 0;
      phiWeight[0] = 1;
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
   void DCHedron::GenerateAngles(uint kind, uint nTrials)
   {
      for (uint i = 0; i < nTrials; ++i)
      {
         data->angles[i] = data->prng.rand(M_PI);
         data->angleEnergy[i] = data->ff.angles.Calc(kind, data->angles[i]);
         data->angleWeights[i] = exp(data->angleEnergy[i] * -data->ff.beta);
      }
   }

   void DCHedron::FreeAngles(uint nTrials)
   {
      for (uint i = 0; i < nBonds; ++i)
      {
         GenerateAngles(angleKinds[i][i], nTrials);
         double stepWeight = std::accumulate(data->angleWeights,
					     data->angleWeights + nTrials, 
					     0.0);
         uint winner = data->prng.PickWeighted(data->angleWeights,
					       nTrials, stepWeight);
         theta[i] = data->angles[winner];
         bendEnergy += data->angleEnergy[winner];
         thetaWeight[i] = stepWeight;
      }
   }

   void DCHedron::PrepareNew()
   {
      bendEnergy = 0;
      FreeAngles(data->nAngleTrials);
      ConstrainedAngles(data->nAngleTrials);
   }

   void DCHedron::PrepareOld()
   {
      bendEnergy = 0;
      FreeAngles(data->nAngleTrials - 1);


   }


   void DCHedron::IncorporateOld(TrialMol& oldMol)
   {
      bendEnergy = 0;
      const Forcefield& ff = data->ff;
      for (uint b = 0; b < nBonds; ++b)
      {
	 
         oldMol.OldThetaAndPhi(bonded[b], focus, theta[b], phi[b]);
         double thetaEnergy = ff.angles.Calc(angleKinds[b][b], theta[b]);
         thetaWeight[b] += exp(-ff.beta * thetaEnergy);
         bendEnergy += thetaEnergy;
	 
	 if (b!=0)
	 {
	    double phiEnergy = 0;
	    phiWeight[b] = 0;
	    for (uint c = 0; c < b; ++c)
	    {
	       double cosTerm = cos(theta[b]) * cos(theta[c]);
	       double sinTerm = sin(theta[b]) * sin(theta[c]);
	       double bfcTheta = acos(sinTerm * cos(phi[b] - phi[c]) + 
				      cosTerm);
	       phiEnergy += ff.angles.Calc(angleKinds[b][c], bfcTheta);
	    }
	    phiWeight[b] = exp(-ff.beta * phiEnergy);
	    bendEnergy += phiEnergy;
	 }
      }
   }


   void DCHedron::ConstrainedAngles(uint nTrials)
   {
      double* angles = data->angles;
      double* energies = data->angleEnergy;
      double* weights = data->angleWeights;
      for (uint b = 1; b < nBonds; ++b)
      {
         //pick "twist" angles 
         for (uint i = 0; i < nTrials; ++i)
	 {
            angles[i] = data->prng.rand(M_PI * 2);
            energies[i] = 0.0;
         }
         //compare to angles determined in previous iterations
         for (uint c = 0; c < b; ++c)
	 {
            double cosTerm = cos(theta[b]) * cos(theta[c]);
            double sinTerm = sin(theta[b]) * sin(theta[c]);
            for (uint i = 0; i < nTrials; ++i)
            {
               double bfcTheta = acos(sinTerm * cos(angles[i] - phi[c]) +
				      cosTerm);
               energies[i] += data->ff.angles.Calc(angleKinds[b][c], bfcTheta);
            }
         }
         //calculate weights from combined energy
         double stepWeight = 0.0;
         for (uint i = 0; i < nTrials; ++i)
	 {
            weights[i] = exp(-1 * data->ff.beta * energies[i]);
            stepWeight += weights[i];
         }
         uint winner = data->prng.PickWeighted(weights, nTrials, stepWeight);
         phi[b] = angles[winner];
         bendEnergy += energies[winner];
         phiWeight[b] = stepWeight;
      }
   }

   //Calculate OldMol Bond Energy &
   //Calculate phi weight for nTrials using actual theta of OldMol
   void DCHedron::ConstrainedAnglesOld(uint nTrials, TrialMol& oldMol)
   {
      IncorporateOld(oldMol);
      
      for (uint b = 1; b < nBonds; ++b) 
      {
	 double stepWeight = 0.0;
	 //pick "twist" angles 
	 for (uint i = 0; i < nTrials; ++i) 
	 {
	    double angles  = data->prng.rand(M_PI * 2);
	    double energies = 0.0;
	    //compare to angles determined in previous iterations
	    for (uint c = 0; c < b; ++c) 
	    {
	       double cosTerm = cos(theta[b]) * cos(theta[c]);
	       double sinTerm = sin(theta[b]) * sin(theta[c]);
	       double bfcTheta = acos(sinTerm * cos(angles - phi[c]) 
				      + cosTerm);
	       energies += data->ff.angles.Calc(angleKinds[b][c], bfcTheta);
	    }
	    
	    //calculate weights from combined energy
	    double weights = exp(-1 * data->ff.beta * energies);
	    stepWeight += weights;
	 }
	 phiWeight[b] += stepWeight;
      }
   }
}

