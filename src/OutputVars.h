#ifndef OUTPUT_VARS_H
#define OUTPUT_VARS_H

#include "BasicTypes.h" //For ulong, uint
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "StaticVals.h"
#include "BoxDimensions.h"
#include "EnergyTypes.h"

class System;
class MoveSettings;
class MoleculeLookup;

class OutputVars
{
public:
   OutputVars(System & sys, StaticVals const& statV);

   ~OutputVars(void);

   void Init(pdb_setup::Atoms const& atoms);
   void InitRef(System & sys, StaticVals const& statV);

   void CalcAndConvert(void);
   uint GetTries(uint sub);
   uint GetAccepted(uint sub);
   double GetAcceptPercent(uint sub);
   double GetScale(uint sub);
   
//private:
   //Intermediate vars.
   uint * numByBox, * numByKindBox;
   double * molFractionByKindBox, * densityByKindBox,
     pressure[BOXES_WITH_U_NB], densityTot[BOX_TOTAL];
   double pTensor[3][3];
   
   uint numKinds;
   //Constants
   double T_in_K;

   //References
   double * volumeRef;
   XYZArray * axisRef;
   double * volInvRef;
   Energy * energyRef, * energyTotRef;
   Virial * virialRef, * virial,  * virialTotRef;
   MoleculeKind * kindsRef;
   MoleculeLookup * molLookupRef;
   
   //Local copy of res names.
   std::vector<std::string> resKindNames;
   double const* movePercRef;
   MoveSettings * moveSetRef;

#if ENSEMBLE == GCMC
   double * chemPot;
#elif ENSEMBLE == GEMC
   double heatOfVap;
#endif
};

#endif /*OUTPUT_VARS_H*/
