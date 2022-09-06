/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/

#include "CheckpointOutput.h"

#include <stdint.h>

#include "GOMC_Config.h"

CheckpointOutput::CheckpointOutput(System &sys, StaticVals const &statV,
                                   Setup &set)
    : moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
      boxDimRef(sys.boxDimRef), molRef(statV.mol), prngRef(sys.prng),
      coordCurrRef(sys.coordinates), trueStepRef(sys.trueStep),
      molSetRef(set.mol), pdbSetupAtomsRef(set.pdb.atoms),
#if GOMC_LIB_MPI
      prngPTRef(*sys.prngParallelTemp),
      enableParallelTemperingBool(sys.ms->parallelTemperingEnabled)
#else
      enableParallelTemperingBool(false)
#endif
{
}

void CheckpointOutput::Init(pdb_setup::Atoms const &atoms,
                            config_setup::Output const &output) {
  enableRestOut = output.restart.settings.enable;
  stepsRestPerOut = output.checkpoint.frequency;
  std::string file = output.statistics.settings.uniqueStr.val + "_restart.chk";
#if GOMC_LIB_MPI
  filename = pathToReplicaOutputDirectory + file;
#else
  filename = file;
#endif
}

void CheckpointOutput::DoOutput(const ulong step) {}

void CheckpointOutput::DoOutputRestart(const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
  std::cout << "Writing checkpoint to file " << filename << " at step "
            << step + 1 << "\n";
  // We want to begin the next simulation on the next step
  // I.e. if we ran 1000 steps, 0-999
  // We want to start on step 1000
  saveCheckpointFile(step + 1, moveSetRef, prngRef, molRef, molLookupRef);
  std::cout << "Checkpoint saved to " << filename << std::endl;
  GOMC_EVENT_STOP(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
}

void CheckpointOutput::saveCheckpointFile(const ulong &step,
                                          MoveSettings &movSetRef, PRNG &prng,
                                          const Molecules &molRef,
                                          MoleculeLookup &molLookRef) {
  std::ofstream ofs(filename);
  if (!ofs.is_open()) {
    fprintf(stderr, "Error writing checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }

  Checkpoint chkObj(
      step, restartFromCheckpoint ? trueStepRef + (step - startStep) : step,
      movSetRef, prng, molRef, molLookRef, molSetRef, pdbSetupAtomsRef);

  cereal::BinaryOutputArchive oa(ofs);
  oa << chkObj;

  oa(cereal::binary_data(molLookRef.molLookup,
                         sizeof(std::uint32_t) * molLookRef.molLookupCount));
  oa(cereal::binary_data(molLookRef.boxAndKindStart,
                         sizeof(std::uint32_t) *
                             molLookRef.boxAndKindStartLength));
  oa(cereal::binary_data(molLookRef.boxAndKindSwappableCounts,
                         sizeof(std::uint32_t) *
                             molLookRef.boxAndKindSwappableLength));
  oa(cereal::binary_data(molLookRef.molIndex,
                         sizeof(std::int32_t) * molLookRef.atomCount));
  oa(cereal::binary_data(molLookRef.atomIndex,
                         sizeof(std::int32_t) * molLookRef.atomCount));
  oa(cereal::binary_data(molLookRef.molKind,
                         sizeof(std::int32_t) * molLookRef.atomCount));
  oa(cereal::binary_data(molLookRef.atomKind,
                         sizeof(std::int32_t) * molLookRef.atomCount));
  oa(cereal::binary_data(molLookRef.atomCharge,
                         sizeof(double) * molLookRef.atomCount));
}
