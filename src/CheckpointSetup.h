/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/

#ifndef CHECKPOINT_SETUP_H
#define CHECKPOINT_SETUP_H

#include <iostream>

#include "Checkpoint.h"
#include "FFSetup.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "MoveSettings.h"
#include "PRNG.h"
#include "Random123Wrapper.h"
#include "VectorLib.h" //for transfer.

class CheckpointSetup {
public:
#if GOMC_LIB_MPI
  CheckpointSetup(ulong &startStep, ulong &trueStep, MoleculeLookup &molLookup,
                  MoveSettings &moveSettings, Molecules &mol, PRNG &prng,
                  Random123Wrapper &r123, Setup &set,
                  const bool &parallelTemperingEnabled, PRNG &prngPT,
                  const std::string &replicaInputDirectoryPath);
#else
  CheckpointSetup(ulong &startStep, ulong &trueStep, MoleculeLookup &molLookup,
                  MoveSettings &moveSettings, Molecules &mol, PRNG &prng,
                  Random123Wrapper &r123, Setup &set);

#endif

  ~CheckpointSetup();

  void loadCheckpointFile();
  void InitOver();

private:
  void SetCheckpointData();

#if GOMC_LIB_MPI
  void SetCheckpointData(const bool &parallelTemperingEnabled, PRNG &prngPT);
#endif

  std::string getFileName();
  void SetStepNumber();
  void SetTrueStepNumber();
  void SetMolecules(Molecules &mols);
  void SetMoveSettings();
  void SetPRNGVariables();
  void SetR123Variables();
  void SetMolecules();
  void SetMoleculeLookup();
  void SetMoleculeSetup();
  void SetPDBSetupAtoms();
#if GOMC_LIB_MPI
  void SetParallelTemperingWasEnabled();
  void SetPRNGVariablesPT(PRNG &prng);
#endif

  void GetOriginalRangeStartStop(uint &_start, uint &stop, const uint m) const;
  void GetRestartRangeStartStop(uint &_start, uint &stop, const uint m) const;

#if GOMC_GTEST || GOMC_GTEST_MPI


#endif

  std::string filename;
  // To avoid repeating Random numbers
  // on the GPU, when InitStep is set to 0
  // we maintain the true step had it not
  // been overwritten by InitStep
  // If init step isn't used
  // trueStep == step
  ulong &startStepRef;
  ulong &trueStepRef;
  MoveSettings &moveSetRef;
  PRNG &prngRef;
  Random123Wrapper &r123Ref;
  Molecules &molRef;
  MoleculeLookup &molLookupRef;
  MolSetup &molSetRef; // 5
  FFSetup &ffSetupRef;
  pdb_setup::Atoms &pdbAtomsRef;

  std::vector<uint> startIdxMolecules;

#if GOMC_LIB_MPI
  bool parallelTemperingWasEnabled;
  const bool &parallelTemperingIsEnabled;
  PRNG &prngPT;
#endif
  Checkpoint chkObj;
  friend class CheckpointOutput;
};

#endif /*CHECKPOINT_SETUP_H*/
