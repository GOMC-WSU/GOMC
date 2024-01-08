/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FREE_ENERGY_OUTPUT_H
#define FREE_ENERGY_OUTPUT_H

#include <fstream>
#include <string>

#include "../lib/BasicTypes.h" //For ulong, uint
#include "../lib/Lambda.h"
#include "../lib/StrLib.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "PDBSetup.h" //For atoms class.
#include "System.h"
#include "UnitConst.h" //For unit conversion factors

struct FreeEnergyOutput : OutputableBase {
  FreeEnergyOutput(OutputVars &v, System &sys);

  ~FreeEnergyOutput();

  virtual void Sample(const ulong step);

  // No additional init.
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  void PrintData(const uint b, const ulong step);
  void CalculateFreeEnergy(const uint b);
  void WriteHeader(void);
  std::string GetString(double a, uint p);

  uint stepsPerSample;
  const CalculateEnergy &calcEn;
  const config_setup::FreeEnergy &freeEnVal;
  const Lambda &lambdaRef;
  Energy dUdL_VDW[BOXES_WITH_U_NB], dUdL_Coulomb[BOXES_WITH_U_NB];
  Energy *energyDiff[BOXES_WITH_U_NB];
  double PV;     // Pressure * Volume
  double Etotal; // Total Energy
  uint lambdaSize, iState;

  std::ofstream outF[BOXES_WITH_U_NB];
  std::string name[BOXES_WITH_U_NB];
#if ENSEMBLE == NPT
  double imposedP; // imposed pressure in NPT
#endif
};

#endif /*FREE_ENERGY_OUTPUT_H*/
