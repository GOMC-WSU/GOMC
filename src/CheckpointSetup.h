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
#include "Random123Wrapper.h"

#include <iostream>
#include "VectorLib.h" //for transfer.

#include "Checkpoint.h"
#include "FFSetup.h"


class CheckpointSetup
{
public:
  CheckpointSetup(ulong & startStep,
                  ulong & trueStep,
                  MoleculeLookup & molLookup, 
                  MoveSettings & moveSettings,
                  Molecules & mol,
                  PRNG & prng,
                  Random123Wrapper & r123,
                  Setup & set);  
#if GOMC_LIB_MPI
                          
  CheckpointSetup(ulong & startStep,
                  ulong & trueStep,
                  MoleculeLookup & molLookup, 
                  MoveSettings & moveSettings,
                  Molecules & mol,
                  PRNG & prng,
                  Random123Wrapper & r123,
                  Setup & set,
                  bool & parallelTemperingEnabled,
                  PRNG & prngPT);
#endif

  ~CheckpointSetup();

  void loadCheckpointFile();
  void InitOver(Molecules & molRef);

private:
  void SetCheckpointData();

  #if GOMC_LIB_MPI
  void SetCheckpointData(bool & parallelTemperingEnabled,
                        PRNG & prngPT);
  #endif

  std::string getFileName();
  void SetStepNumber();
  void SetTrueStepNumber();
  void SetMolecules(Molecules& mols);
  void SetMoleculeKindDictionary(Molecules& mols);
  void SetMoveSettings();
  void SetPRNGVariables();
  void SetR123Variables();
  void SetMolecules();
  void SetMoleculeKindDictionary();
  void SetMoleculeIndices();
  void SetMoleculeSetup();
  void SetPDBSetupAtoms();
#if GOMC_LIB_MPI  
  void SetParallelTemperingWasEnabled();
  void SetPRNGVariablesPT(PRNG & prng);
#endif

#if GOMC_GTEST

#endif

  std::string filename;
  // To avoid repeating Random numbers
  // on the GPU, when InitStep is set to 0
  // we maintain the true step had it not
  // been overwritten by InitStep
  // If init step isn't used
  // trueStep == step
  ulong & startStepRef;
  ulong & trueStepRef;
  MoveSettings & moveSetRef;
  PRNG & prngRef;
  Random123Wrapper & r123Ref;
  Molecules & molRef;
  MoleculeLookup & molLookupRef;
  MolSetup & molSetRef;        //5
  FFSetup & ffSetupRef;
  pdb_setup::Atoms & pdbAtomsRef;

#if GOMC_LIB_MPI
  bool parallelTemperingWasEnabled;
  bool & parallelTemperingIsEnabled;
  PRNG & prngPT;
#endif
  Checkpoint chkObj;  
  friend class CheckpointOutput;
};

#endif
