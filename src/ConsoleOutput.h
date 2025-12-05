/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H

#include "BasicTypes.h" //For uint
#include "MoleculeKind.h"
#include "Molecules.h"
#include "MoveConst.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "PDBSetup.h"
#include "StaticVals.h"

class System;
namespace config_setup {
struct Output;
}
class SystemPotential;
class Energy;
class Virial;
class MoveSettings;
class MoleculeLookup;

struct ConsoleOutput : OutputableBase {
public:
  ConsoleOutput(OutputVars &v) { this->var = &v; }

  // Console Output does not need to sample, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output) {
    enableOut = output.console.enable;
    stepsPerOut = output.console.frequency;
    enableEnergy = output.statistics.vars.energy.fluct;
    enablePressure = output.statistics.vars.pressure.fluct;
    enableSurfTension = output.statistics.vars.surfaceTension.fluct;

#ifdef VARIABLE_VOLUME
    enableVolume = output.statistics.vars.volume.fluct;
#else
    enableVolume = false;
#endif

#ifdef VARIABLE_PARTICLE_NUMBER
    enableMol = output.statistics.vars.molNum.fluct;
#else
    enableMol = false;
#endif
    enableDens = output.statistics.vars.density.fluct;
    if (enableVolume || enablePressure || enableMol || enableDens ||
        enableSurfTension) {
      enableStat = true;
    }
    WriteConsoleHeaders = true;
    DoOutput(0);
  }
  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  const static int elementWidth = 21;
  bool enableEnergy, enablePressure, enableDens, enableVolume, enableMol;
  bool enableSurfTension, enableStat;
  bool WriteConsoleHeaders;
  void PrintMove(const uint box, const ulong step) const;
  void PrintMoveStat(const uint box, const ulong step) const;
  void PrintStatistic(const uint box, const ulong step) const;
  void PrintPressureTensor(const uint box, const ulong step) const;
  void PrintEnergy(const uint box, Energy const &en, const ulong step) const;
  void PrintEnergyTitle();
  void PrintStatisticTitle();
  void PrintMoveTitle();
  void printElement(const double t, const int width, uint precision = 8) const;
  void printElement(const uint t, const int width) const;
  void printElement(const std::string t, const int width) const;

  template <typename T>
  void printElementStep(const T t, const ulong step, const int width) const;
};

#endif /*CONSOLE_OUTPUT_H*/
