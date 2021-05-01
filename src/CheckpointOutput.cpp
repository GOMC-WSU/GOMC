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

namespace
{
union dbl_output_union {
  char bin_value[8];
  double dbl_value;
};

union uint32_output_union {
  char bin_value[4];
  uint32_t uint_value;
};

union uint64_output_union {
  char bin_value[8];
  uint64_t uint_value;
};

union int8_input_union {
  char bin_value[1];
  int8_t int_value;
};
}

CheckpointOutput::CheckpointOutput(System & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), prngRef(sys.prng),
  coordCurrRef(sys.coordinates),
#if GOMC_LIB_MPI
  prngPTRef(*sys.prngParallelTemp),
  enableParallelTempering(sys.ms->parallelTemperingEnabled)
#else
  enableParallelTempering(false)
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
  openOutputFile();
  printGOMCVersion();
  printStepNumber(step);
  printRandomNumbers();
  printMoleculeLookupData();
  printMoveSettingsData();
  printMoleculesData();
#if GOMC_LIB_MPI
  printParallelTemperingBoolean();
  if(enableParallelTempering)
    printRandomNumbersParallelTempering();
#endif
  std::cout << "Checkpoint saved to " << filename << std::endl;
  GOMC_EVENT_STOP(1, GomcProfileEvent::CHECKPOINT_OUTPUT);
}

void CheckpointOutput::setGOMCVersion()
{
  sprintf(gomc_version, "%d.%02d", GOMC_VERSION_MAJOR, GOMC_VERSION_MINOR % 100);
}

void CheckpointOutput::printGOMCVersion()
{
  setGOMCVersion();
  fprintf(outputFile, "%c%s%c", '$', gomc_version, '$');
}

void CheckpointOutput::printParallelTemperingBoolean()
{
  int8_t s = (int8_t) enableParallelTempering;
  write_uint8_binary(s);
}

void CheckpointOutput::printStepNumber(const ulong step)
{
  write_uint64_binary(step+1);
}

void CheckpointOutput::printRandomNumbers()
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  const int N = 624;
  // The MT::save function also appends the "left" variable,
  // so need to allocate one more array element
  uint32_t* saveArray = new uint32_t[N + 1];
  prngRef.GetGenerator()->save(saveArray);
  for(int i = 0; i < N; i++) {
    write_uint32_binary(saveArray[i]);
  }

  // Save the location of pointer in state
  uint32_t location = prngRef.GetGenerator()->pNext -
                      prngRef.GetGenerator()->state;
  write_uint32_binary(location);

  // save the "left" value so we can restore it later
  write_uint32_binary(prngRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  write_uint32_binary(prngRef.GetGenerator()->seedValue);

  delete[] saveArray;
}

#if GOMC_LIB_MPI
void CheckpointOutput::printRandomNumbersParallelTempering()
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  const int N = 624;
  uint32_t* saveArray = new uint32_t[N];
  prngPTRef.GetGenerator()->save(saveArray);
  for(int i = 0; i < N; i++) {
    write_uint32_binary(saveArray[i]);
  }

  // Save the location of pointer in state
  uint32_t location = prngPTRef.GetGenerator()->pNext -
                      prngPTRef.GetGenerator()->state;
  write_uint32_binary(location);

  // save the "left" value so we can restore it later
  write_uint32_binary(prngPTRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  write_uint32_binary(prngPTRef.GetGenerator()->seedValue);
}
#endif

void CheckpointOutput::printMoleculeLookupData()
{
  // print the size of molLookup array
  write_uint32_binary(molLookupRef.molLookupCount);
  // print the molLookup array itself
  for(int i = 0; i < (int) molLookupRef.molLookupCount; i++) {
    write_uint32_binary(molLookupRef.molLookup[i]);
  }

  // print the size of boxAndKindStart array
  write_uint32_binary(molLookupRef.boxAndKindStartCount);
  // print the BoxAndKindStart array
  for(int i = 0; i < (int) molLookupRef.boxAndKindStartCount; i++) {
    write_uint32_binary(molLookupRef.boxAndKindStart[i]);
  }

  // print numKinds
  write_uint32_binary(molLookupRef.numKinds);
  //print the size of fixedMolecule array
  write_uint32_binary((uint)molLookupRef.fixedMolecule.size());
  //print the fixedMolecule array itself
  for(int i = 0; i < (int) molLookupRef.fixedMolecule.size(); i++) {
    write_uint32_binary(molLookupRef.fixedMolecule[i]);
  }
}

void CheckpointOutput::printMoleculesData()
{
  // print the start of each molecule
  // there is an extra one at the end which store the total count
  // that is used for to calculate the length of molecule
  // so the length of last molecule can be calculated using
  // start[molIndex+1]-start[molIndex]
  write_uint32_binary(molRef.count);
  for(int i = 0; i < (int)molRef.count+1; i++) {
    write_uint32_binary(molRef.start[i]);
  }

  // print the start of each kind
  write_uint32_binary(molRef.kIndexCount);
  for(int i = 0; i < (int)molRef.kIndexCount; i++) {
    write_uint32_binary(molRef.kIndex[i]);
  }
}

void CheckpointOutput::printMoveSettingsData()
{
  printVector3DDouble(moveSetRef.scale);
  printVector3DDouble(moveSetRef.acceptPercent);
  printVector3DUint(moveSetRef.accepted);
  printVector3DUint(moveSetRef.tries);
  printVector3DUint(moveSetRef.tempAccepted);
  printVector3DUint(moveSetRef.tempTries);
  printVector2DUint(moveSetRef.mp_tries);
  printVector2DUint(moveSetRef.mp_accepted);
  printVector1DDouble(moveSetRef.mp_t_max);
  printVector1DDouble(moveSetRef.mp_r_max);
}

void CheckpointOutput::printVector3DDouble(const std::vector< std::vector< std::vector<double> > > &data)
{
  // print size of tempTries
  ulong size_x = data.size();
  ulong size_y = data[0].size();
  ulong size_z = data[0][0].size();
  write_uint64_binary(size_x);
  write_uint64_binary(size_y);
  write_uint64_binary(size_z);

  // print tempTries array
  for(int i = 0; i < (int) size_x; i++) {
    for(int j = 0; j < (int) size_y; j++) {
      for(int k = 0; k < (int) size_z; k++) {
        write_double_binary(data[i][j][k]);
      }
    }
  }
}

void CheckpointOutput::printVector3DUint(const std::vector< std::vector< std::vector<uint> > > &data)
{
  // print size of tempTries
  ulong size_x = data.size();
  ulong size_y = data[0].size();
  ulong size_z = data[0][0].size();
  write_uint64_binary(size_x);
  write_uint64_binary(size_y);
  write_uint64_binary(size_z);

  // print tempTries array
  for(int i = 0; i < (int) size_x; i++) {
    for(int j = 0; j < (int) size_y; j++) {
      for(int k = 0; k < (int) size_z; k++) {
        write_uint32_binary(data[i][j][k]);
      }
    }
  }
}

void CheckpointOutput::printVector2DUint(const std::vector< std::vector< uint > > &data)
{
  // print size of array
  ulong size_x = data.size();
  ulong size_y = data[0].size();
  write_uint64_binary(size_x);
  write_uint64_binary(size_y);

  // print array itself
  for(int i = 0; i < (int) size_x; i++) {
    for(int j = 0; j < (int) size_y; j++) {
      write_uint32_binary(data[i][j]);
    }
  }
}

void CheckpointOutput::printVector1DDouble(const std::vector< double > &data)
{
  // print size of array
  ulong size_x = data.size();
  write_uint64_binary(size_x);

  // print array itself
  for(int i = 0; i < (int) size_x; i++) {
    write_double_binary(data[i]);
  }
}

void CheckpointOutput::openOutputFile()
{
  outputFile = fopen(filename.c_str(), "wb");
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
}

void CheckpointOutput::write_double_binary(double data)
{
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  dbl_output_union temp;
  temp.dbl_value = data;
  fprintf(outputFile, "%c%c%c%c%c%c%c%c",
          temp.bin_value[0],
          temp.bin_value[1],
          temp.bin_value[2],
          temp.bin_value[3],
          temp.bin_value[4],
          temp.bin_value[5],
          temp.bin_value[6],
          temp.bin_value[7]);
  fflush(outputFile);
}

void CheckpointOutput::write_uint8_binary(int8_t data)
{
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  int8_input_union temp;
  temp.int_value = data;
  fprintf(outputFile, "%c",
          temp.bin_value[0]);
  fflush(outputFile);
}

void CheckpointOutput::write_uint64_binary(uint64_t data)
{
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint64_output_union temp;
  temp.uint_value = htof64(data);
  fprintf(outputFile, "%c%c%c%c%c%c%c%c",
          temp.bin_value[0],
          temp.bin_value[1],
          temp.bin_value[2],
          temp.bin_value[3],
          temp.bin_value[4],
          temp.bin_value[5],
          temp.bin_value[6],
          temp.bin_value[7]);
  fflush(outputFile);
}

void CheckpointOutput::write_uint32_binary(uint32_t data)
{
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint32_output_union temp;
  temp.uint_value = htof32(data);
  fprintf(outputFile, "%c%c%c%c",
          temp.bin_value[0],
          temp.bin_value[1],
          temp.bin_value[2],
          temp.bin_value[3]);
  fflush(outputFile);
}
