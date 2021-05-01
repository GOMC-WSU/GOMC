/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef OUTPUT_ABSTRACTS_H
#define OUTPUT_ABSTRACTS_H

#include <string>

#include "BasicTypes.h" //For ulong

#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "StaticVals.h"
#include "ConfigSetup.h" //For enables, etc.
#include "PDBSetup.h" //For atoms class
#include "GOMCEventsProfile.h" // for profiling

#include "GOMC_Config.h"    //For MPI 
#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

class OutputVars;
class System;

////////////////////////////
///  WRITE
//

class OutputableBase
{
public:
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output) = 0;

  virtual void DoOutput(const ulong step) = 0;

  virtual void DoOutputRestart(const ulong step) = 0;

  virtual void Sample(const ulong step) = 0;

  virtual void Output(const ulong step) final
  {
    /* merged psf only prints on first step */
    if (!enableOut && !enableRestOut && !forceOutput) {
      return;
    } else {
      Sample(step);
    }

    /* We will output either when the step number is every stepsPerOut
       Or recalculate trajectory is enabled (forceOutput) */
    if ((onlyPrintOnFirstStep && firstPrint) || (enableOut && ((step + 1) % stepsPerOut == 0)) || forceOutput) {
    /*
    I want to always generate a merged psf file, even on recalcTraj where startStep != 0
    //if ((onlyPrintOnFirstStep && step == 0) || (enableOut && ((step + 1) % stepsPerOut == 0)) || forceOutput) {
    */
      DoOutput(step);
      firstPrint = false;
    }

    /* We will output if the step number is every stepsRestPerOut */
    /* (forceOutput && onlyPrintOnFirstStep) - only ever true for PSF files, which I want to see to test collation */
    if ((enableRestOut && ((step + 1) % stepsRestPerOut == 0)) || (forceOutput && onlyPrintOnFirstStep)) {
      DoOutputRestart(step);
    }
  }

  void Init(pdb_setup::Atoms const& atoms,
            config_setup::Output const& output,
            const ulong tillEquil,
            const ulong totSteps)
  {
    Init(tillEquil, totSteps, output.statistics.settings.uniqueStr.val);
    Init(atoms, output);
  }

  void Init(const ulong tillEquil, const ulong totSteps,
            std::string const& uniqueForFileIO)
  {
#if GOMC_LIB_MPI
    // For some reason, whether or not split passes by reference or value changes per OS.
    // Hence the need for a disposable copy of the string
    std::string copyOfUniqueForFileIO = uniqueForFileIO.c_str();
#ifdef WIN32
    std::vector<std::string> tokens = OutputableBase::split(copyOfUniqueForFileIO, std::string(2, OS_SEP));
#else
    std::vector<std::string> tokens = OutputableBase::split(copyOfUniqueForFileIO, std::string(1, OS_SEP));
#endif
    std::stringstream replicaDirectory;
    // Loop terminates before last entry to remove the value that
    // The user gets passes as "OutputName", which we have prefixed w a path
    // In the ParallelTemperingPreproccessing in Main.cpp.
    // This leaves only the path, which we store as "pathToReplicaOutputDirectory"
    for(int i = 0; i < tokens.size() - 1; ++i) {
      replicaDirectory << tokens[i] << OS_SEP;
    }
    pathToReplicaOutputDirectory = replicaDirectory.str();

    uniqueName = tokens[tokens.size() - 1];
#else
    uniqueName = uniqueForFileIO;
#endif
    stepsTillEquil = tillEquil;
    totSimSteps = totSteps;

    // We will use forceOutput for recalculate trajectory
    // If we are not running any simulation then the step will stay 0
    // But we still need to output console for each frame
    // So by setting this value to true we will force all the outputs to
    // output even though the step is zero
    if(totSteps == 0)
      forceOutput = true;
    else
      forceOutput = false;
  }

#if GOMC_LIB_MPI
  std::vector<std::string> split(std::string str, std::string sep)
  {
    char* cstr = const_cast<char*>(str.c_str());
    char* current;
    std::vector<std::string> arr;
    current = strtok(cstr, sep.c_str());
    while(current != NULL) {
      arr.push_back(current);
      current = strtok(NULL, sep.c_str());
    }
    return arr;
  }
#endif

//private:
  std::string uniqueName;
#if GOMC_LIB_MPI
  std::string pathToReplicaOutputDirectory;
#endif
  ulong stepsPerOut = 0, stepsRestPerOut = 0, stepsTillEquil = 0, totSimSteps = 0;
  bool enableOut = false, enableRestOut = false, firstPrint = true, forceOutput = false, onlyPrintOnFirstStep = false;

  //Contains references to various objects.
  OutputVars * var;
};

#endif /*OUTPUT_ABSTRACTS_H*/
