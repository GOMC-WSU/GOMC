/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FREEENERGY_OUTPUT_H
#define FREEENERGY_OUTPUT_H

#include <string>
#include <fstream>

#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "../lib/BasicTypes.h" //For ulong, uint
#include "../lib/StrLib.h"
#include "../lib/Lambda.h"
#include "System.h"
#include "PDBSetup.h" //For atoms class.
#include "EnergyTypes.h"
#include "CalculateEnergy.h"
#include "UnitConst.h" //For unit conversion factors

struct FreeEnergyOutput : OutputableBase {

  FreeEnergyOutput(OutputVars & v, System & sys);

  ~FreeEnergyOutput();

  virtual void Sample(const ulong step);

  //No additional init.
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void DoOutput(const ulong step);

private:
  void PrintData(const uint b, const uint step);
  void CalculateFreeEnergy(const uint b);
  void WriteHeader(void);
  std::string GetString(double a, uint p);

  uint stepsPerSample;
  const CalculateEnergy& calcEn;
  const config_setup::FreeEnergy&  freeEnVal;
  const Lambda& lambdaRef;
  Energy dUdL_VDW[BOXES_WITH_U_NB], dUdL_Coulomb[BOXES_WITH_U_NB];
  Energy *energyDiff[BOXES_WITH_U_NB];
  double PV; // Pressure * Volume
  double Etotal; //Total Energy
  uint lambdaSize, iState;

  std::ofstream outF[BOXES_WITH_U_NB];
  std::string name[BOXES_WITH_U_NB];

};

#endif /*HIST_OUTPUT_H*/
