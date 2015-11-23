#include "OutputVars.h"
#include "UnitConst.h" //For unit conversion factors
#include "System.h"
#include "MoveSettings.h"
#include "MoleculeLookup.h"

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

OutputVars::OutputVars(System & sys, StaticVals const& statV): molFractionByKindBox(NULL)
{ InitRef(sys, statV); }

void OutputVars::InitRef(System & sys, StaticVals const& statV)
{
   T_in_K = statV.forcefield.T_in_K;
   volumeRef = sys.boxDimRef.volume;
   axisRef = &sys.boxDimRef.axis;
   volInvRef = sys.boxDimRef.volInv;
   energyTotRef = & sys.potential.totalEnergy;
   virialTotRef = & sys.potential.totalVirial;
   energyRef = sys.potential.boxEnergy;
   virialRef = sys.potential.boxVirial;
   kindsRef = statV.mol.kinds;
   molLookupRef = & sys.molLookupRef;
   moveSetRef = & sys.moveSettings;
   movePercRef = statV.movePerc;
   
   virial = new Virial[BOXES_WITH_U_NB];
}

uint OutputVars::GetTries(uint sub)
{
   return (sub < mv::SCALEABLE ?
           moveSetRef->tries[sub] + moveSetRef->tempTries[sub] :
           moveSetRef->tries[sub]);
}

uint OutputVars::GetAccepted(uint sub)
{
   return (sub < mv::SCALEABLE ?
           moveSetRef->accepted[sub] + moveSetRef->tempAccepted[sub] :
           moveSetRef->accepted[sub]);
}

double OutputVars::GetScale(uint sub) { return moveSetRef->scale[sub]; }

double OutputVars::GetAcceptPercent(uint sub)
{ return (double)(GetAccepted(sub))/(double)(GetTries(sub))*100.0; }

void OutputVars::Init(pdb_setup::Atoms const& atoms)
{
   //Init vals.
   numKinds = molLookupRef->GetNumKind();
   resKindNames = atoms.resKindNames;

   //Allocate arrays,
   uint kTot = BOX_TOTAL * numKinds;
   numByBox = new uint [BOX_TOTAL];
   numByKindBox = new uint [kTot];
   densityByKindBox = new double [kTot];
   if (numKinds > 1)
      molFractionByKindBox = new double [kTot];
}

OutputVars::~OutputVars(void)
{
   if ( numByBox != NULL )
      delete[] numByBox;
   if ( numByKindBox != NULL )
      delete[] numByKindBox;
   if ( molFractionByKindBox != NULL )
      delete[] molFractionByKindBox;
   if ( densityByKindBox != NULL )
      delete[] densityByKindBox;
}

void OutputVars::CalcAndConvert(void)
{
   double rawPressure[BOXES_WITH_U_NB];

#if ENSEMBLE == GEMC
   //Which box is the liquid in gibbs ensemble
   uint liqBox = mv::BOX0, vapBox = mv::BOX1;
#endif

   molLookupRef->TotalAndDensity(numByBox,  numByKindBox, molFractionByKindBox,
				 densityByKindBox, volInvRef);

#if ENSEMBLE == GEMC
   //Determine which box is liquid for purposes of heat of vap.
   if (densityByKindBox[numKinds] > densityByKindBox[0])
   {
      vapBox = mv::BOX0;
      liqBox = mv::BOX1;
   }
#endif

   for (uint b = 0; b < BOXES_WITH_U_NB; b++)
   {
      //Account for dimensionality of virial (raw "virial" is actually a
      //multiple of the true virial, based on the dimensions stress is exerted
      //in)
      virial[b] = virialRef[b];
      virial[b] /= unit::DIMENSIONALITY;
      virial[b] /= volumeRef[b];
      rawPressure[b] = 0.0;
      for (uint k = 0; k < numKinds; k++)
      {
	 double * density = &densityByKindBox[k+numKinds*b];

	 // Instead of converting to mass first
	 // (an alternate route to calculate the ideal gas pressure)
	 // the form of the Boltzmann constant that kcal/mol/K is used
	 // such that a single conversion factor can be applied to both
	 // the ideal and virial components of the pressure.
	 rawPressure[b] += *density;
      }
      // Finish ideal component
      rawPressure[b] *= T_in_K;
      // Add the virial component
      rawPressure[b] -= virial[b].total;

      // Convert to desired units
      // ( starting: K * molecule / Angstrom^3 )
      pressure[b] = rawPressure[b];
      pressure[b] *= unit::K_MOLECULE_PER_A3_TO_BAR;
   }
   for (uint b = 0; b < BOX_TOTAL; b++)
   {
      for (uint k = 0; k < numKinds; k++)
      {
	 double * density = &densityByKindBox[k+numKinds*b];

	 // Convert density to g/ml (which is equivalent to g/cm3)
	 // To get kg/m3, multiply output densities by 1000.
	 *density *= unit::MOLECULES_PER_A3_TO_MOL_PER_CM3 * 
	    kindsRef[k].molMass;
      }
   }
#if ENSEMBLE == GEMC
   // delta Hv = (Uv-Ul) + P(Vv-Vl)
   heatOfVap = energyRef[vapBox].total/numByBox[vapBox] - 
      energyRef[liqBox].total/numByBox[liqBox] +
      rawPressure[vapBox]*(volumeRef[vapBox]/numByBox[vapBox] - 
			   volumeRef[liqBox]/numByBox[liqBox]);
   heatOfVap *= unit::K_TO_KJ_PER_MOL;
#endif

}
