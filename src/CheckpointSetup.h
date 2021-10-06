/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef CHECKPOINT_SETUP_H
#define CHECKPOINT_SETUP__H

#include "MoveSettings.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "PRNG.h"

#include <iostream>
#include "VectorLib.h" //for transfer.

#include "Checkpoint.h"


class CheckpointSetup
{
public:
  CheckpointSetup(MoleculeLookup & molLookup, 
                                MoveSettings & moveSettings,
                                Molecules & mol,
                                PRNG & prng,
                                Setup const& set);
  ~CheckpointSetup();

  void loadCheckpointFile(ulong & startStep);
  void InitOver(Molecules & molRef);

private:
  void SetCheckpointData   (ulong & startStep,
                            MoveSettings & movSetRef,
                            PRNG & prng,
                            Molecules & molRef,
                            MoleculeLookup & molLookRef);

  void SetCheckpointData   (ulong & startStep);

  #if GOMC_LIB_MPI
  void SetCheckpointData   (ulong & startStep,
                            MoveSettings & movSetRef,
                            PRNG & prng,
                            Molecules & molRef,
                            MoleculeLookup & molLookRef)
                            bool & parallelTemperingEnabled,
                            PRNG & prngPT);
  #endif

  std::string getFileName();
  void SetStepNumber(ulong & startStep);
  void SetTrueStepNumber(ulong & trueStep);
  void SetMoveSettings(MoveSettings & movSetRef);
  void SetPRNGVariables(PRNG & prng);
  void SetMolecules(Molecules& mols);
  void SetMoleculeKindDictionary(Molecules& mols);
  void SetMoleculeIndices(MoleculeLookup& molLookup);
  void SetMoveSettings();
  void SetPRNGVariables();
  void SetMolecules();
  void SetMoleculeKindDictionary();
  void SetMoleculeIndices();
#if GOMC_LIB_MPI  
  void SetParallelTemperingWasEnabled();
  void SetPRNGVariablesPT(PRNG & prng);
#endif

  std::string filename;
  // To avoid repeating Random numbers
  // on the GPU, when InitStep is set to 0
  // we maintain the true step had it not
  // been overwritten by InitStep
  // If init step isn't used
  // trueStep == step
  ulong  trueStepNumber;
  MoveSettings & moveSetRef;
  PRNG & prngRef;
  Molecules & molRef;
  MoleculeLookup & molLookupRef;
#if GOMC_LIB_MPI
  bool parallelTemperingWasEnabled;
  bool & parallelTemperingIsEnabled;
  PRNG & prngPT;
#endif
  Checkpoint chkObj;  
  friend class CheckpointOutput;
};

#endif
