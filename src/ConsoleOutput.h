#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H

#include "BasicTypes.h" //For uint
#include "OutputAbstracts.h"
#include "Molecules.h"
#include "MoleculeKind.h"
#include "StaticVals.h"
#include "PDBSetup.h"
#include "MoveConst.h"
#include "OutputVars.h"

class System;
namespace config_setup
{
struct Output;
}
class SystemPotential;
class Energy;
class Virial;
class MoveSettings;
class MoleculeLookup;

struct ConsoleOutput : OutputableBase
{
public:
  ConsoleOutput(OutputVars & v)
  {
    this->var = &v;
  }

  //Console Output does not need to sample, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output)
  {
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
	enableSurfTension)
    {
      enableStat = true;
    }
    DoOutput(0);
  }
  virtual void DoOutput(const ulong step);

private:
  const static int elementWidth = 16;
  bool enableEnergy, enablePressure, enableDens, enableVolume, enableMol;
  bool enableSurfTension, enableStat;
  void PrintMove(const uint box, const ulong step) const;
  void PrintMoveStat(const uint box, const ulong step) const;
  void PrintStatistic(const uint box) const;
  void PrintPressureTensor(const uint box) const;
  void PrintEnergy(const uint box, Energy const& en, Virial const& vir) const;
  void PrintEnergyTitle(const uint box);
  void PrintStatisticTitle(const uint box);
  void PrintMoveTitle(const uint box);
  template <typename T> void printElement ( const T t, const int width) const;
};

#endif /*CONSOLE_OUTPUT_H*/
