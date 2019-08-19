/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef OUTPUT_VARS_H
#define OUTPUT_VARS_H

#include "BasicTypes.h" //For ulong, uint
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "StaticVals.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "EnergyTypes.h"
#include "CalculateEnergy.h"

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

  void CalcAndConvert(ulong step);
  bool Performed(uint moveKind);
  uint GetTries(uint box, uint sub);
  uint GetAccepted(uint box, uint sub);
  double GetAcceptPercent(uint box, uint sub);
  double GetScale(uint box, uint sub);

//private:
  //Intermediate vars.
  uint * numByBox, * numByKindBox;
  double * molFractionByKindBox, * densityByKindBox,
         pressure[BOXES_WITH_U_NB], densityTot[BOX_TOTAL];
  double pressureTens[BOXES_WITH_U_NB][3][3];
  double surfaceTens[BOXES_WITH_U_NB];
  ulong pCalcFreq;
  bool pressureCalc;

  uint numKinds;
  //Constants
  double T_in_K;

  //References
  double * volumeRef;
  XYZArray axisRef;
  double * volInvRef;
  Energy * energyRef, * energyTotRef;
  Virial * virialRef, * virial,  * virialTotRef;
  MoleculeKind * kindsRef;
  MoleculeLookup * molLookupRef;
  CalculateEnergy& calc;

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
