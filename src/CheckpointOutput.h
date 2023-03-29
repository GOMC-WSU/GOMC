/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/

#ifndef CHECKPOINT_OUTPUT_H
#define CHECKPOINT_OUTPUT_H

#include <iostream>

#include "Checkpoint.h"
#include "Coordinates.h"
#include "GOMC_Config.h"
#include "MoveBase.h"
#include "MoveSettings.h"
#include "OutputAbstracts.h"

class CheckpointOutput : public OutputableBase {
public:
  CheckpointOutput(System &sys, StaticVals const &statV, Setup &set);

  ~CheckpointOutput() {}

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);
  virtual void Sample(const ulong step) {}

private:
  void saveCheckpointFile(const ulong &startStep, MoveSettings &movSetRef,
                          PRNG &prng, const Molecules &molRef,
                          MoleculeLookup &molLookRef);
#if GOMC_LIB_MPI
  void saveCheckpointFile(const ulong &startStep, MoveSettings &movSetRef,
                          PRNG &prng, const Molecules &molRef,
                          MoleculeLookup &molLookRef,
                          bool &parallelTemperingEnabled, PRNG &prngPT);
#endif

  // To avoid repeating Random numbers
  // on the GPU, when InitStep is set to 0
  // we maintain the true step had it not
  // been overwritten by InitStep
  // If init step isn't used
  // trueStep == step
  ulong &trueStepRef;
  MoveSettings &moveSetRef;
  MoleculeLookup &molLookupRef;
  BoxDimensions &boxDimRef;
  Molecules const &molRef;
  PRNG &prngRef;
  Coordinates &coordCurrRef;
  MolSetup &molSetRef; // 5
  pdb_setup::Atoms &pdbSetupAtomsRef;
#if GOMC_LIB_MPI
  PRNG &prngPTRef;
#endif

  bool enableParallelTemperingBool;
  std::string filename;
  ulong stepsPerCheckpoint;
};

#endif
