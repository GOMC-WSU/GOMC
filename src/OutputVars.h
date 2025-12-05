/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef OUTPUT_VARS_H
#define OUTPUT_VARS_H

#include "BasicTypes.h" //For ulong, uint
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "StaticVals.h"

class System;
class MoveSettings;
class MoleculeLookup;

class OutputVars {
public:
  OutputVars(System &sys, StaticVals const &statV,
             const std::vector<std::string> &molKindNames);

  ~OutputVars(void);

  void Init();
  void InitRef(System &sys, StaticVals const &statV);

  void CalcAndConvert(ulong step);
  bool Performed(uint moveKind);
  uint GetTries(uint box, uint sub);
  uint GetAccepted(uint box, uint sub);
  double GetAcceptPercent(uint box, uint sub);
  double GetScale(uint box, uint sub);

  // private:
  // Intermediate vars.
  uint *numByBox, *numByKindBox;
  double *molFractionByKindBox, *densityByKindBox, pressure[BOXES_WITH_U_NB],
      densityTot[BOX_TOTAL], compressibility[BOXES_WITH_U_NB],
      enthalpy[BOXES_WITH_U_NB];
  double pressureTens[BOXES_WITH_U_NB][3][3];
  double surfaceTens[BOXES_WITH_U_NB];
  ulong pCalcFreq;
  bool pressureCalc;

  uint numKinds;
  // Constants
  const double &T_in_K;

  // References
  double *volumeRef;
  XYZArray *axisRef;
  double *volInvRef;
  Energy *energyRef, *energyTotRef;
  Virial *virialRef, *virial, *virialTotRef;
  MoleculeKind *kindsRef;
  MoleculeLookup *molLookupRef;
  CalculateEnergy &calc;

  // Local copy of res names.
  std::vector<std::string> molKindNames;
  double const *movePercRef;
  MoveSettings *moveSetRef;

#if ENSEMBLE == GCMC
  double *chemPot;
#elif ENSEMBLE == GEMC
  // Which box is the liquid in gibbs ensemble
  uint liqBox, vapBox;
  double heatOfVap;
  double heatOfVap_energy_term_box[BOX_TOTAL];
  double heatOfVap_density_term_box[BOX_TOTAL];
#endif
};

#endif /*OUTPUT_VARS_H*/
