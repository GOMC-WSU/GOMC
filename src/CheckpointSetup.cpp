/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
  real dbl_value;
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

  for(int b = 0; b < totalBoxes; b++) {
    axis[b].resize(3);
    cosAngle[b].resize(3);
    axis[b][0] = readRealIn8Chars();
    axis[b][1] = readRealIn8Chars();
    axis[b][2] = readRealIn8Chars();
    cosAngle[b][0] = readRealIn8Chars();
    cosAngle[b][1] = readRealIn8Chars();
    cosAngle[b][2] = readRealIn8Chars();
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
  for(int i = 0; i < coordLength; i++) {
    XYZ temp;
    temp.x = readRealIn8Chars();
    temp.y = readRealIn8Chars();
    temp.z = readRealIn8Chars();
    coords.Set(i, temp);
  }
}

void CheckpointSetup::readMoleculeLookupData()
{
  // read the size of molLookup array
  molLookupVec.resize(readUintIn8Chars());

  // read the molLookup array itself
  for(int i = 0; i < molLookupVec.size(); i++) {
    molLookupVec[i] = readUintIn8Chars();
  }

  // read the size of boxAndKindStart array
  boxAndKindStartVec.resize(readUintIn8Chars());
  // read the BoxAndKindStart array
  for(int i = 0; i < boxAndKindStartVec.size(); i++) {
    boxAndKindStartVec[i] = readUintIn8Chars();
  }

  // read numKinds
  numKinds = readUintIn8Chars();

  //read the size of fixedAtom array
  fixedAtomVec.resize(readUintIn8Chars());
  //read the fixedAtom array itself
  for(int i = 0; i < fixedAtomVec.size(); i++) {
    fixedAtomVec[i] = readUintIn8Chars();
  }
}

void CheckpointSetup::readMoveSettingsData()
{
  uint size_x, size_y, size_z;

  // read size of scale
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // read scale array
  scaleVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    scaleVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      scaleVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        scaleVec[i][j][k] = readRealIn8Chars();
      }
    }
  }

  // read size of acceptPercent
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // read acceptPercent array
  acceptPercentVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    acceptPercentVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      acceptPercentVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        acceptPercentVec[i][j][k] = readRealIn8Chars();
      }
    }
  }

  // read size of accepted
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // read accepted array
  acceptedVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    acceptedVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      acceptedVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        acceptedVec[i][j][k] = readUintIn8Chars();
      }
    }
  }

  // print size of tries
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // print tries array
  triesVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    triesVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      triesVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        triesVec[i][j][k] = readUintIn8Chars();
      }
    }
  }

  // print size of tempAccepted
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // print tempAccepted array
  tempAcceptedVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    tempAcceptedVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      tempAcceptedVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        tempAcceptedVec[i][j][k] = readUintIn8Chars();
      }
    }
  }

  // print size of tempTries
  size_x = readUintIn8Chars();
  size_y = readUintIn8Chars();
  size_z = readUintIn8Chars();

  // print tempTries array
  tempTriesVec.resize(size_x);
  for(int i = 0; i < size_x; i++) {
    tempTriesVec[i].resize(size_y);
    for(int j = 0; j < size_y; j++) {
      tempTriesVec[i][j].resize(size_z);
      for(int k = 0; k < size_z; k++) {
        tempTriesVec[i][j][k] = readUintIn8Chars();
      }
    }
  }
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

real CheckpointSetup::readRealIn8Chars()
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

void CheckpointSetup::SetStepNumber(ulong & startStep)
{
  startStep = stepNumber;
}

void CheckpointSetup::SetBoxDimensions(BoxDimensions & boxDimRef)
{
  for(int b = 0; b < totalBoxes; b++) {
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
  for(int i = 0; i < this->molLookupVec.size(); i++) {
    molLookupRef.molLookup[i] = this->molLookupVec[i];
  }
  for(int i = 0; i < this->boxAndKindStartVec.size(); i++) {
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
}
