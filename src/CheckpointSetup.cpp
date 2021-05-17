/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointSetup.h"
#include "MoleculeLookup.h"
#include "System.h"

#include "Endian.h"

namespace
{
union dbl_input_union {
  char bin_value[8];
  double dbl_value;
};

union uint64_input_union {
  char bin_value[8];
  uint64_t uint_value;
};

union uint32_input_union {
  char bin_value[4];
  uint32_t uint_value;
};

union int8_input_union {
  char bin_value[1];
  int8_t int_value;
};
}

CheckpointSetup::CheckpointSetup(System & sys, StaticVals const& statV,
                                 Setup const& set) :
  parallelTemperingWasEnabled(false)
{
  std::string file = set.config.in.files.checkpoint.name[0];
#if GOMC_LIB_MPI
  filename = sys.ms->replicaInputDirectoryPath + file;
#else
  filename = file;
#endif
  inputFile = NULL;
  saveArray = NULL;
  #if GOMC_LIB_MPI
  saveArrayPT = NULL;
  #endif
}

CheckpointSetup::CheckpointSetup(std::string file) :
  parallelTemperingWasEnabled(false)
{
  filename = file;
  inputFile = NULL;
  saveArray = NULL;
  #if GOMC_LIB_MPI
  saveArrayPT = NULL;
  #endif
}

void CheckpointSetup::ReadAll()
{
  openInputFile();
  readGOMCVersion();
  readStepNumber();
  readRandomNumbers();
  readMoveSettingsData();
  //readMoleculeLookupData();
  readMoleculesData();
  readOriginalMoleculeIndices();
#if GOMC_LIB_MPI
  readParallelTemperingBoolean();
  if(parallelTemperingWasEnabled)
    readRandomNumbersParallelTempering();
#endif
  std::cout << "Checkpoint loaded from " << filename << std::endl;
}

bool CheckpointSetup::isLegacy()
{
  char first_symbol = ' ';
  int ret = fscanf(inputFile, "%c", &first_symbol);
  if(ret != 1) {
    std::cerr << "Could not read checkpoint file!\n";
    exit(EXIT_FAILURE);
  }
  return first_symbol != '$';
}

void CheckpointSetup::readGOMCVersion()
{
  if(isLegacy()) {
    // Return cursor to beginning of the file
    fseek(inputFile, 0, SEEK_SET);
    sprintf(gomc_version, "0.00");
  } else {
    int result = fscanf(inputFile, "%[^$]", gomc_version);
    // Move the cursor past the $ sign
    fseek(inputFile, 1, SEEK_CUR);
  }
}

void CheckpointSetup::readParallelTemperingBoolean()
{
  parallelTemperingWasEnabled = read_uint8_binary();
}

void CheckpointSetup::readStepNumber()
{
  stepNumber = read_uint64_binary();
}

void CheckpointSetup::readRandomNumbers()
{
  // First let's read the state array
  // the length of the array is 624
  const int N = 624;
  if(saveArray != NULL) {
    delete[] saveArray;
    saveArray = NULL;
  }
  saveArray = new uint32_t[N];
  for(int i = 0; i < N; i++) {
    saveArray[i] = read_uint32_binary();
  }

  // Read the location of pointer in state
  seedLocation = read_uint32_binary();

  // Read the "left" value so we can restore
  seedLeft = read_uint32_binary();

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedValue = read_uint32_binary();
}

#if GOMC_LIB_MPI

void CheckpointSetup::readRandomNumbersParallelTempering()
{
  // First let's read the state array
  // the length of the array is 624
  const int N = 624;
  if(saveArrayPT != NULL) {
    delete[] saveArray;
    saveArrayPT = NULL;
  }
  saveArrayPT = new uint32_t[N];
  for(int i = 0; i < N; i++) {
    saveArrayPT[i] = read_uint32_binary();
  }

  // Read the location of pointer in state
  seedLocationPT = read_uint32_binary();

  // Read the "left" value so we can restore
  seedLeftPT = read_uint32_binary();

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedValuePT = read_uint32_binary();
}

#endif

void CheckpointSetup::readMoveSettingsData()
{
  readVector3DDouble(scaleVec);
  readVector3DDouble(acceptPercentVec);
  readVector3DUint(acceptedVec);
  readVector3DUint(triesVec);
  readVector3DUint(tempAcceptedVec);
  readVector3DUint(tempTriesVec);
  readVector2DUint(mp_triesVec);
  readVector2DUint(mp_acceptedVec);
  readVector1DDouble(mp_t_maxVec);
  readVector1DDouble(mp_r_maxVec);
}

void CheckpointSetup::readMoleculesData()
{
  // read the size of start array
  uint startCount = read_uint32_binary() + 1;
  molecules_originalStartVec.resize(startCount);
  for(int i = 0; i < (int)startCount; i++) {
    molecules_originalStartVec[i] = read_uint32_binary();
  }

  // read the kIndex array
  molecules_originalKIndexVec.resize(read_uint32_binary());
  for(int i = 0; i < (int) molecules_originalKIndexVec.size(); i++) {
    molecules_originalKIndexVec[i] = read_uint32_binary();
  }
}

void CheckpointSetup::readOriginalMoleculeIndices(){
 readVector1DUint(originalMoleculeIndicesVec);
 readVector1DUint(permutedMoleculeIndicesVec);
}

void CheckpointSetup::openInputFile()
{
  inputFile = fopen(filename.c_str(), "rb");
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint input file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
}

void CheckpointSetup::closeInputFile()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Checkpoint file was not open!\n");
    exit(EXIT_FAILURE);
  }
  fclose(inputFile);
}

void CheckpointSetup::SetStepNumber(ulong & startStep)
{
  startStep = stepNumber;
}

void CheckpointSetup::SetPRNGVariables(PRNG & prng)
{
  prng.GetGenerator()->load(saveArray);
  prng.GetGenerator()->pNext = prng.GetGenerator()->state + seedLocation;
  prng.GetGenerator()->left = seedLeft;
  prng.GetGenerator()->seedValue = seedValue;
}

bool CheckpointSetup::CheckIfParallelTemperingWasEnabled()
{
  return (bool)parallelTemperingWasEnabled;
}


#if GOMC_LIB_MPI
void CheckpointSetup::SetPRNGVariablesPT(PRNG & prng)
{
  prng.GetGenerator()->load(saveArrayPT);
  prng.GetGenerator()->pNext = prng.GetGenerator()->state + seedLocationPT;
  prng.GetGenerator()->left = seedLeftPT;
  prng.GetGenerator()->seedValue = seedValuePT;
}
#endif

void CheckpointSetup::SetMoveSettings(MoveSettings & moveSettings)
{
  moveSettings.scale = this->scaleVec;
  moveSettings.acceptPercent = this->acceptPercentVec;
  moveSettings.accepted = this->acceptedVec;
  moveSettings.tries = this->triesVec;
  moveSettings.tempAccepted = this->tempAcceptedVec;
  moveSettings.tempTries = this->tempTriesVec;
  moveSettings.mp_tries = this->mp_triesVec;
  moveSettings.mp_accepted = this->mp_acceptedVec;
  moveSettings.mp_t_max = this->mp_t_maxVec;
  moveSettings.mp_r_max = this->mp_r_maxVec;
}

void CheckpointSetup::SetMolecules(Molecules& mols)
{
  for(int i = 0; i < (int)this->molecules_originalStartVec.size(); i++) {
    mols.originalStart[i] = molecules_originalStartVec[i];
  }
  for(int i = 0; i < (int)this->molecules_originalKIndexVec.size(); i++) {
    mols.originalKIndex[i] = molecules_originalKIndexVec[i];
  }
}

/* Index magic. Don't touch */
void CheckpointSetup::SetOriginalMoleculeIndices(MoleculeLookup & molLookupRef)
{
  /* Original Mol Indices are for constant trajectory output from start to finish of a single run*/
  molLookupRef.originalMoleculeIndices = vect::transfer<uint32_t>(this->originalMoleculeIndicesVec);
  /* Permuted Mol Indices are for following single molecules as molLookup permutes the indices and continuing the next run*/
  molLookupRef.permutedMoleculeIndices = vect::transfer<uint32_t>(this->permutedMoleculeIndicesVec);
}

void
CheckpointSetup::readVector3DDouble(std::vector<std::vector<std::vector<double> > > &data)
{
  // read size of data
  ulong size_x = read_uint64_binary();
  ulong size_y = read_uint64_binary();
  ulong size_z = read_uint64_binary();

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j].resize(size_z);
      for(int k = 0; k < (int) size_z; k++) {
        data[i][j][k] = read_double_binary();
      }
    }
  }
}

void CheckpointSetup::readVector3DUint(std::vector<std::vector<std::vector<uint> > > &data)
{
  // read size of data
  ulong size_x = read_uint64_binary();
  ulong size_y = read_uint64_binary();
  ulong size_z = read_uint64_binary();

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j].resize(size_z);
      for(int k = 0; k < (int) size_z; k++) {
        data[i][j][k] = read_uint32_binary();
      }
    }
  }
}

void CheckpointSetup::readVector2DUint(std::vector<std::vector<uint> > &data)
{
  // read size of data
  ulong size_x = read_uint64_binary();
  ulong size_y = read_uint64_binary();

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i].resize(size_y);
    for(int j = 0; j < (int) size_y; j++) {
      data[i][j] = read_uint32_binary();
    }
  }
}

void CheckpointSetup::readVector1DDouble(std::vector<double> &data)
{
// read size of data
  ulong size_x = read_uint64_binary();

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i] = read_double_binary();
  }
}

void CheckpointSetup::readVector1DUint(std::vector< uint > &data)
{
// read size of data
  ulong size_x = read_uint64_binary();

  // read array
  data.resize(size_x);
  for(int i = 0; i < (int) size_x; i++) {
    data[i] = read_uint32_binary();
  }
}
