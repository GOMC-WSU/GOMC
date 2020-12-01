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
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), coordCurrRef(sys.coordinates),
  prngRef(sys.prng)
{
  std::string file = set.config.in.files.checkpoint.name[0];
#if GOMC_LIB_MPI
  filename = sys.ms->replicaInputDirectoryPath + file;
#else
  filename = file;
#endif
  inputFile = NULL;
  saveArray = NULL;
}

void CheckpointSetup::ReadAll()
{
  openInputFile();
  readGOMCVersion();
  readStepNumber();
  readRandomNumbers();
  readMoleculeLookupData();
  readMoveSettingsData();
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
  return first_symbol != '$';
}

void CheckpointSetup::readGOMCVersion()
{
  if(isLegacy()) {
    // Return cursor to beginning of the file
    fseek(inputFile, 0, SEEK_SET);
    sprintf(gomc_version, "0.00\0");
  } else {
    fscanf(inputFile, "%[^$]", gomc_version);
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

void CheckpointSetup::readBoxDimensionsData()
{
  // read the number of boxes
  totalBoxes = read_uint32_binary();
  axis.resize(totalBoxes);
  cosAngle.resize(totalBoxes);

  for(int b = 0; b < (int) totalBoxes; b++) {
    axis[b].resize(3);
    cosAngle[b].resize(3);
    axis[b][0] = read_double_binary();
    axis[b][1] = read_double_binary();
    axis[b][2] = read_double_binary();
    cosAngle[b][0] = read_double_binary();
    cosAngle[b][1] = read_double_binary();
    cosAngle[b][2] = read_double_binary();
  }
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

void CheckpointSetup::readCoordinates()
{
  // first let's read the count
  coordLength = read_uint32_binary();

  // now let's read the coordinates one by one
  coords.Init(coordLength);
  for(int i = 0; i < (int) coordLength; i++) {
    XYZ temp;
    temp.x = read_double_binary();
    temp.y = read_double_binary();
    temp.z = read_double_binary();
    coords.Set(i, temp);
  }
}

void CheckpointSetup::readMoleculeLookupData()
{
  // read the size of molLookup array
  molLookupVec.resize(read_uint32_binary());
  // read the molLookup array itself
  for(int i = 0; i < (int) molLookupVec.size(); i++) {
    molLookupVec[i] = read_uint32_binary();
  }

  // read the size of boxAndKindStart array
  boxAndKindStartVec.resize(read_uint32_binary());
  // read the BoxAndKindStart array
  for(int i = 0; i < (int) boxAndKindStartVec.size(); i++) {
    boxAndKindStartVec[i] = read_uint32_binary();
  }

  // read numKinds
  numKinds = read_uint32_binary();
  //read the size of fixedAtom array
  fixedAtomVec.resize(read_uint32_binary());
  //read the fixedAtom array itself
  for(int i = 0; i < (int) fixedAtomVec.size(); i++) {
    fixedAtomVec[i] = read_uint32_binary();
  }
}

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

void CheckpointSetup::openInputFile()
{
  inputFile = fopen(filename.c_str(), "rb");
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint input file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
}

double CheckpointSetup::read_double_binary()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  dbl_input_union temp;
  int ret = fscanf(inputFile, "%c%c%c%c%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3],
                   &temp.bin_value[4],
                   &temp.bin_value[5],
                   &temp.bin_value[6],
                   &temp.bin_value[7]);
  if(ret != 8) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  return temp.dbl_value;
}

int8_t CheckpointSetup::read_uint8_binary()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  int8_input_union temp;
  int ret = fscanf(inputFile, "%c",
                   &temp.bin_value[0]);
  if(ret != 1) {
    // We could add this back if we REQUIRE the PT flag as output, but
    // this would break all previous checkpoint files generated by legacy code.
    //std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    //exit(EXIT_FAILURE);
    // If we ran out of binary, then return 0, which evaluates to false
    // when casted to boolean.  This way legacy checkpoints can be used
    // as input to MPI builds.
    return 0;
  }
  return temp.int_value;
}

uint32_t CheckpointSetup::read_uint32_binary()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint32_input_union temp;
  int ret = fscanf(inputFile, "%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3]);
  if(ret != 4) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  // Fix endianness, implementation in lib/Endian.h
  temp.uint_value = ftoh32(temp.uint_value);
  return temp.uint_value;
}

uint64_t CheckpointSetup::read_uint64_binary()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint64_input_union temp;
  int ret = fscanf(inputFile, "%c%c%c%c%c%c%c%c",
                   &temp.bin_value[0],
                   &temp.bin_value[1],
                   &temp.bin_value[2],
                   &temp.bin_value[3],
                   &temp.bin_value[4],
                   &temp.bin_value[5],
                   &temp.bin_value[6],
                   &temp.bin_value[7]);
  if(ret != 8) {
    std::cerr << "CheckpointSetup couldn't read required data from binary!\n";
    exit(EXIT_FAILURE);
  }
  // Fix endianness, implementation in lib/Endian.h
  temp.uint_value = ftoh64(temp.uint_value);
  return temp.uint_value;
}

void CheckpointSetup::SetStepNumber(ulong & startStep)
{
  startStep = stepNumber;
}

void CheckpointSetup::SetBoxDimensions(BoxDimensions & boxDimRef)
{
  for(int b = 0; b < (int) totalBoxes; b++) {
    boxDimRef.axis.Set(b, axis[b][0], axis[b][1], axis[b][2]);
    boxDimRef.cosAngle[b][0] = this->cosAngle[b][0];
    boxDimRef.cosAngle[b][1] = this->cosAngle[b][1];
    boxDimRef.cosAngle[b][2] = this->cosAngle[b][2];
    boxDimRef.volume[b] = axis[b][0] * axis[b][1] * axis[b][2];
    boxDimRef.volInv[b] = 1.0 / boxDimRef.volume[b];
  }
  boxDimRef.axis.CopyRange(boxDimRef.halfAx, 0, 0, BOX_TOTAL);
  boxDimRef.halfAx.ScaleRange(0, BOX_TOTAL, 0.5);
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

void CheckpointSetup::SetCoordinates(Coordinates & coordinates)
{
  coords.CopyRange(coordinates, 0, 0, coordLength);
}

void CheckpointSetup::SetMoleculeLookup(MoleculeLookup & molLookupRef)
{
  if(molLookupRef.molLookupCount != this->molLookupVec.size()) {
    std::cerr << "ERROR: Restarting from checkpoint...\n"
              << "molLookup size does not match with restart file\n";
    exit(EXIT_FAILURE);
  }
  for(int i = 0; i < (int) this->molLookupVec.size(); i++) {
    molLookupRef.molLookup[i] = this->molLookupVec[i];
  }
  for(int i = 0; i < (int) this->boxAndKindStartVec.size(); i++) {
    molLookupRef.boxAndKindStart[i] = this->boxAndKindStartVec[i];
  }
  molLookupRef.numKinds = this->numKinds;
}

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
