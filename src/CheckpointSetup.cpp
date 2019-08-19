/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointSetup.h"
#include "MoleculeLookup.h"
#include "System.h"

namespace
{
  union dbl_input_union {
    char bin_value[8];
    double dbl_value;
  };

  union uint32_input_union {
    char bin_value[8];
    uint32_t uint_value;
  };
}

CheckpointSetup::CheckpointSetup(System & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), prngRef(sys.prng),
  coordCurrRef(sys.coordinates), filename("checkpoint.dat")
{
  inputFile = NULL;
  saveArray = NULL;
}

void CheckpointSetup::ReadAll()
{
  openInputFile();
  readStepNumber();
  readBoxDimensionsData();
  readRandomNumbers();
  readCoordinates();
  readMoleculeLookupData();
  readMoveSettingsData();
  std::cout << "Checkpoint loaded from " << filename << std::endl;
}

void CheckpointSetup::readStepNumber()
{
  stepNumber = readUintIn8Chars();
}

void CheckpointSetup::readBoxDimensionsData()
{
  // read the number of boxes
  totalBoxes = readUintIn8Chars();
  axis.resize(totalBoxes);
  cosAngle.resize(totalBoxes);

  for(int b=0; b<totalBoxes; b++) {
    axis[b].resize(3);
    cosAngle[b].resize(3);
    axis[b][0] = readDoubleIn8Chars();
    axis[b][1] = readDoubleIn8Chars();
    axis[b][2] = readDoubleIn8Chars();
    cosAngle[b][0] = readDoubleIn8Chars();
    cosAngle[b][1] = readDoubleIn8Chars();
    cosAngle[b][2] = readDoubleIn8Chars();
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
  for(int i=0; i<N; i++) {
    saveArray[i] = readUintIn8Chars();
  }

  // Read the location of pointer in state
  seedLocation = readUintIn8Chars();

  // Read the "left" value so we can restore
  seedLeft = readUintIn8Chars();

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedValue = readUintIn8Chars();
}

void CheckpointSetup::readCoordinates()
{
  // first let's read the count
  coordLength = readUintIn8Chars();

  // now let's read the coordinates one by one
  coords.Init(coordLength);
  for(int i=0; i<coordLength; i++) {
    XYZ temp;
    temp.x = readDoubleIn8Chars();
    temp.y = readDoubleIn8Chars();
    temp.z = readDoubleIn8Chars();
    coords.Set(i, temp);
  }
}

void CheckpointSetup::readMoleculeLookupData()
{
  // read the size of molLookup array
  molLookupVec.resize(readUintIn8Chars());

  // read the molLookup array itself
  for(int i=0; i<molLookupVec.size(); i++) {
    molLookupVec[i] = readUintIn8Chars();
  }

  // read the size of boxAndKindStart array
  boxAndKindStartVec.resize(readUintIn8Chars());
  // read the BoxAndKindStart array
  for(int i=0; i<boxAndKindStartVec.size(); i++) {
    boxAndKindStartVec[i] = readUintIn8Chars();
  }

  // read numKinds
  numKinds = readUintIn8Chars();

  //read the size of fixedAtom array
  fixedAtomVec.resize(readUintIn8Chars());
  //read the fixedAtom array itself
  for(int i=0; i<fixedAtomVec.size(); i++) {
    fixedAtomVec[i] = readUintIn8Chars();
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

double CheckpointSetup::readDoubleIn8Chars()
{
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  dbl_input_union temp;
  fscanf(inputFile, "%c%c%c%c%c%c%c%c",
         &temp.bin_value[0],
         &temp.bin_value[1],
         &temp.bin_value[2],
         &temp.bin_value[3],
         &temp.bin_value[4],
         &temp.bin_value[5],
         &temp.bin_value[6],
         &temp.bin_value[7]);
  return temp.dbl_value;
}

uint32_t CheckpointSetup::readUintIn8Chars()
{
  uint32_t data;
  if(inputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint32_input_union temp;
  fscanf(inputFile, "%c%c%c%c%c%c%c%c",
         &temp.bin_value[0],
         &temp.bin_value[1],
         &temp.bin_value[2],
         &temp.bin_value[3],
         &temp.bin_value[4],
         &temp.bin_value[5],
         &temp.bin_value[6],
         &temp.bin_value[7]);
  return temp.uint_value;
}

void CheckpointSetup::SetStepNumber(ulong & startStep) {
  startStep = stepNumber;
}

void CheckpointSetup::SetBoxDimensions(BoxDimensions & boxDimRef) {
  for(int b=0; b<totalBoxes; b++) {
    boxDimRef.axis[b][0] = axis[b][0];
    boxDimRef.axis[b][1] = axis[b][1];
    boxDimRef.axis[b][2] = axis[b][2];
    boxDimRef.cosAngle[b][0] = this->cosAngle[b][0];
    boxDimRef.cosAngle[b][1] = this->cosAngle[b][1];
    boxDimRef.cosAngle[b][2] = this->cosAngle[b][2];
    boxDimRef.volume[b] = axis[b][0] * axis[b][1] * axis[b][2];
    boxDimRef.volInv[b] = 1.0 / boxDimRef.volume[b];
  }

  for (uint b = 0; b < BOX_TOTAL; b++) {
    boxDimRef.halfAx[b][0] = boxDimRef.axis[b][0] * 0.5;
    boxDimRef.halfAx[b][1] = boxDimRef.axis[b][1] * 0.5;
    boxDimRef.halfAx[b][2] = boxDimRef.axis[b][2] * 0.5;
  }
}

void CheckpointSetup::SetPRNGVariables(PRNG & prng) {
  prng.GetGenerator()->load(saveArray);
  prng.GetGenerator()->pNext = prng.GetGenerator()->state + seedLocation;
  prng.GetGenerator()->left = seedLeft;
  prng.GetGenerator()->seedValue = seedValue;
}

void CheckpointSetup::SetCoordinates(Coordinates & coordinates) {
  coords.CopyRange(coordinates, 0, 0, coordLength);
}

void CheckpointSetup::SetMoleculeLookup(MoleculeLookup & molLookupRef) {
  if(molLookupRef.molLookupCount != this->molLookupVec.size()) {
    std::cerr << "ERROR: Restarting from checkpoint...\n"
      << "molLookup size does not match with restart file\n";
    exit(EXIT_FAILURE);
  }
  for(int i=0; i<this->molLookupVec.size(); i++) {
    molLookupRef.molLookup[i] = this->molLookupVec[i];
  }
  for(int i=0; i<this->boxAndKindStartVec.size(); i++) {
    molLookupRef.boxAndKindStart[i] = this->boxAndKindStartVec[i];
  }
  molLookupRef.numKinds = this->numKinds;
}

void CheckpointSetup::SetMoveSettings(MoveSettings & moveSettings) {
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
CheckpointSetup::readVector3DDouble(vector<vector<vector<double> > > &data) {
  // read size of data
  ulong size_x = readUintIn8Chars();
  ulong size_y = readUintIn8Chars();
  ulong size_z = readUintIn8Chars();

  // read array
  data.resize(size_x);
  for(int i=0; i<size_x; i++) {
    data[i].resize(size_y);
    for(int j=0; j<size_y; j++) {
      data[i][j].resize(size_z);
      for(int k=0; k<size_z; k++) {
        data[i][j][k] = readDoubleIn8Chars();
      }
    }
  }
}

void CheckpointSetup::readVector3DUint(vector<vector<vector<uint> > > &data) {
  // read size of data
  ulong size_x = readUintIn8Chars();
  ulong size_y = readUintIn8Chars();
  ulong size_z = readUintIn8Chars();

  // read array
  data.resize(size_x);
  for(int i=0; i<size_x; i++) {
    data[i].resize(size_y);
    for(int j=0; j<size_y; j++) {
      data[i][j].resize(size_z);
      for(int k=0; k<size_z; k++) {
        data[i][j][k] = readUintIn8Chars();
      }
    }
  }
}

void CheckpointSetup::readVector2DUint(vector<vector<uint> > &data) {
  // read size of data
  ulong size_x = readUintIn8Chars();
  ulong size_y = readUintIn8Chars();

  // read array
  data.resize(size_x);
  for(int i=0; i<size_x; i++) {
    data[i].resize(size_y);
    for(int j=0; j<size_y; j++) {
      data[i][j] = readUintIn8Chars();
    }
  }
}

void CheckpointSetup::readVector1DDouble(vector<double> &data) {
// read size of data
  ulong size_x = readUintIn8Chars();

  // read array
  data.resize(size_x);
  for(int i=0; i<size_x; i++) {
    data[i] = readDoubleIn8Chars();
  }
}
