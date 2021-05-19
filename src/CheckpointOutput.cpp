/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointOutput.h"
#include "MoleculeLookup.h"
#include "System.h"
#include "GOMC_Config.h"

#include "Endian.h"

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
  outputFile = NULL;
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
  //openOutputFile();
  setGOMCVersion();
  setStepNumber(step);
  setRandomNumbers();
  /* For consistent trajectory ordering */
  setSortedMoleculeIndices();
#if GOMC_LIB_MPI
  setParallelTemperingBoolean();
  if(enableParallelTempering)
    setRandomNumbersParallelTempering();
#endif
  // create and open a character archive for output
  openOutputFile(filename);
  *oa << *this;
  closeOutputFile();
  std::cout << "Checkpoint saved to " << filename << std::endl;
  GOMC_EVENT_STOP(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
}

void CheckpointOutput::openOutputFile(std::string filenameArg)
{
  ofs = new std::ofstream(filenameArg);
  oa = new boost::archive::text_oarchive(*ofs);
  if(!ofs->is_open()) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filenameArg.c_str());
    exit(EXIT_FAILURE);
  }
}

void CheckpointOutput::closeOutputFile()
{
  if(!ofs->is_open()) {
    fprintf(stderr, "Checkpoint file was not open!\n");
    exit(EXIT_FAILURE);
  }
  ofs->close();
}

void CheckpointOutput::setGOMCVersion()
{
  sprintf(gomc_version, "%d.%02d", GOMC_VERSION_MAJOR, GOMC_VERSION_MINOR % 100);
}

void CheckpointOutput::setParallelTemperingBoolean()
{
  enableParallelTempering = (int8_t) enableParallelTemperingBool;
}

void CheckpointOutput::setStepNumber(const ulong stepArg)
{
   step = (uint64_t) stepArg;
}

void CheckpointOutput::setRandomNumbers()
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  // The MT::save function also appends the "left" variable,
  // so need to allocate one more array element
  prngRef.GetGenerator()->save(saveArray);
  // Save the location of pointer in state
  location = prngRef.GetGenerator()->pNext - prngRef.GetGenerator()->state;

  // save the "left" value so we can restore it later
  left = (prngRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seed = prngRef.GetGenerator()->seedValue;
}

#if GOMC_LIB_MPI
void CheckpointOutput::setRandomNumbersParallelTempering()
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  // The MT::save function also appends the "left" variable,
  // so need to allocate one more array element
  prngPTRef.GetGenerator()->save(saveArrayPT);
  // Save the location of pointer in state
  locationPT = prngPTRef.GetGenerator()->pNext - prngPTRef.GetGenerator()->state;

  // save the "left" value so we can restore it later
  leftPT = (prngPTRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedPT = prngPTRef.GetGenerator()->seedValue;
}
#endif

/* After the first run, the molecules are sorted, so we need to use the same sorting process
   seen below, to reinitialize the originalMolInds every checkpoint */
void CheckpointOutput::setSortedMoleculeIndices(){
  uint molCounter = 0, b, k, kI, countByKind, molI;
  if (!molLookupRef.restartFromCheckpoint){
    for (b = 0; b < BOX_TOTAL; ++b) {
      for (k = 0; k < molLookupRef.numKinds; ++k) {
        countByKind = molLookupRef.NumKindInBox(k, b);
        for (kI = 0; kI < countByKind; ++kI) {
          molI = molLookupRef.GetMolNum(kI, k, b);
          originalMoleculeIndicesCheck[molCounter] = molI;
          ++molCounter;
        }
      }
    }
    for (uint molI = 0; molI < molLookupRef.molLookupCount; ++molI){
      molLookupRef.permutedMoleculeIndices[molI] = molLookupRef.originalMoleculeIndices[molI];
    }
  } else {
    for (b = 0; b < BOX_TOTAL; ++b) {
      for (k = 0; k < molLookupRef.numKinds; ++k) {
        countByKind = molLookupRef.NumKindInBox(k, b);
        for (kI = 0; kI < countByKind; ++kI) {
          molI = molLookupRef.GetSortedMolNum(kI, k, b);
          originalMoleculeIndicesCheck[molCounter] = molLookupRef.permutedMoleculeIndices[molI];
          ++molCounter;
        }
      }
    }
  }
}
