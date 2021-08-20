/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointOutput.h"
#include "GOMC_Config.h"


CheckpointOutput::CheckpointOutput(System & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), prngRef(sys.prng),
  coordCurrRef(sys.coordinates),
#if GOMC_LIB_MPI
  prngPTRef(*sys.prngParallelTemp),
  enableParallelTemperingBool(sys.ms->parallelTemperingEnabled)
#else
  enableParallelTemperingBool(false)
#endif
{
}

void CheckpointOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output)
{
  enableRestOut = output.restart.settings.enable;
  stepsRestPerOut = output.restart.settings.frequency;
  std::string file = output.statistics.settings.uniqueStr.val + "_restart.chk";
#if GOMC_LIB_MPI
  filename = pathToReplicaOutputDirectory + file;
#else
  filename = file;
#endif
}

void CheckpointOutput::DoOutput(const ulong step){}

void CheckpointOutput::DoOutputRestart(const ulong step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
  std::cout << "Writing checkpoint to file " << filename << " at step " << step+1 << "\n";
  saveCheckpointFile(step, moveSetRef, prngRef, molRef, molLookupRef);
  std::cout << "Checkpoint saved to " << filename << std::endl;
  GOMC_EVENT_STOP(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
}

void CheckpointOutput::saveCheckpointFile(const ulong & startStep,
                                          MoveSettings & movSetRef,
                                          PRNG & prng,
                                          const Molecules & molRef,
                                          MoleculeLookup & molLookRef){
  std::ofstream ofs(filename);
  if (!ofs.is_open()){
    fprintf(stderr, "Error writing checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }

  Checkpoint chkObj(startStep,
                    movSetRef,
                    prng,
                    molRef,
                    molLookRef);

  #if GOMC_BOOST_LIB
    boost::archive::text_oarchive oa(ofs);
    oa << chkObj;
  #else
    cereal::BinaryOutputArchive oa(ofs);
    oa << chkObj;
  #endif
}




