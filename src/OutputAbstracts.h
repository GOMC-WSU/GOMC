/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef OUTPUT_ABSTRACTS_H
#define OUTPUT_ABSTRACTS_H

#include <string>

#include "BasicTypes.h"           //For ulong
#include "ConfigSetup.h"          //For enables, etc.
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "GOMCEventsProfile.h"    // for profiling
#include "GOMC_Config.h"          //For MPI
#include "PDBSetup.h"             //For atoms class
#include "StaticVals.h"
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

class OutputableBase {
public:
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output) = 0;

  virtual void DoOutput(const ulong step) = 0;

  virtual void DoOutputRestart(const ulong step) = 0;

  virtual void Sample(const ulong step) = 0;

  virtual void Output(const ulong step) final {
    if ((!enableOut && !enableRestOut && !forceOutput)) {
      return;
    } else {
      Sample(step);
    }

    /* We will output either when the step number is every stepsPerOut
       Or recalculate trajectory is enabled (forceOutput) */
    /* printOnFirstStep -- only true for PSFOutput, WolfCalibration  */
    if ((printOnFirstStep && step == startStep) ||
        (enableOut && ((step + 1) % stepsPerOut == 0) || forceOutput)) {
      DoOutput(step);
      firstPrint = false;
    }

    /* We will output if the step number is every stepsRestPerOut */
    if ((enableRestOut && ((step + 1) % stepsRestPerOut == 0)) || forceOutput) {
      DoOutputRestart(step);
    }
  }

  void Init(pdb_setup::Atoms const &atoms, config_setup::Input const &input,
            config_setup::Output const &output,
            config_setup::SystemVals const &sys, const ulong startStep,
            const ulong tillEquil, const ulong totSteps) {
    initStepRead = sys.step.initStepRead;
    this->startStep = startStep;
    restartFromCheckpoint = input.restart.restartFromCheckpoint;
    Init(tillEquil, totSteps, output.statistics.settings.uniqueStr.val);
    Init(atoms, output);
  }

  void Init(const ulong tillEquil, const ulong totSteps,
            std::string const &uniqueForFileIO) {
#if GOMC_LIB_MPI
    // For some reason, whether or not split passes by reference or value
    // changes per OS. Hence the need for a disposable copy of the string
    std::string copyOfUniqueForFileIO = uniqueForFileIO.c_str();
#ifdef WIN32
    std::vector<std::string> tokens =
        OutputableBase::split(copyOfUniqueForFileIO, std::string(2, OS_SEP));
#else
    std::vector<std::string> tokens =
        OutputableBase::split(copyOfUniqueForFileIO, std::string(1, OS_SEP));
#endif
    std::stringstream replicaDirectory;
    // Loop terminates before last entry to remove the value that
    // The user gets passes as "OutputName", which we have prefixed w a path
    // In the ParallelTemperingPreproccessing in Main.cpp.
    // This leaves only the path, which we store as
    // "pathToReplicaOutputDirectory"
    for (int i = 0; i < tokens.size() - 1; ++i) {
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
    if (totSteps == 0)
      forceOutput = true;
    else
      forceOutput = false;
  }

#if GOMC_LIB_MPI
  std::vector<std::string> split(std::string str, std::string sep) {
    char *cstr = const_cast<char *>(str.c_str());
    char *current;
    std::vector<std::string> arr;
    current = strtok(cstr, sep.c_str());
    while (current != NULL) {
      arr.push_back(current);
      current = strtok(NULL, sep.c_str());
    }
    return arr;
  }
#endif

  // private:
  std::string uniqueName;
#if GOMC_LIB_MPI
  std::string pathToReplicaOutputDirectory;
#endif
  ulong stepsPerOut = 0, stepsRestPerOut = 0, stepsTillEquil = 0,
        totSimSteps = 0, startStep = 0;
  bool enableOut = false, enableRestOut = false, firstPrint = true,
       forceOutput = false, printOnFirstStep = false;

  // Contains references to various objects.
  OutputVars *var;
  bool restartFromCheckpoint, initStepRead;
};

#endif /*OUTPUT_ABSTRACTS_H*/
