/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointOutput.h"
#include "MoleculeLookup.h"
#include "System.h"

namespace
{
  union dbl_output_union {
    char bin_value[8];
    double dbl_value;
  };

  union uint32_output_union {
    char bin_value[8];
    uint32_t uint_value;
  };
}

CheckpointOutput::CheckpointOutput(System & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), prngRef(sys.prng),
  coordCurrRef(sys.coordinates), filename("checkpoint.dat")
{
  outputFile = NULL;
}

void CheckpointOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output)
{
  enableOutCheckpoint = output.checkpoint.enable;
  stepsPerCheckpoint = output.checkpoint.frequency;
}

void CheckpointOutput::DoOutput(const ulong step)
{
  if(enableOutCheckpoint) {
    openOutputFile();
    printStepNumber(step);
    printBoxDimensionsData();
    printRandomNumbers();
    printCoordinates();
    printMoleculeLookupData();
    printMoveSettingsData();
    std::cout << "Checkpoint saved to " << filename << std::endl;
  }
}

void CheckpointOutput::printStepNumber(const ulong step)
{
  uint32_t s = (uint32_t) step + 1;
  outputUintIn8Chars(s);
}

void CheckpointOutput::printBoxDimensionsData()
{
  // print the number of boxes
  uint32_t totalBoxes = BOX_TOTAL;
  outputUintIn8Chars(totalBoxes);
  for(int b=0; b<totalBoxes; b++) {
    XYZ axis = boxDimRef.axis.Get(b);
    outputDoubleIn8Chars(axis.x);
    outputDoubleIn8Chars(axis.y);
    outputDoubleIn8Chars(axis.z);
    outputDoubleIn8Chars(boxDimRef.cosAngle[b][0]);
    outputDoubleIn8Chars(boxDimRef.cosAngle[b][1]);
    outputDoubleIn8Chars(boxDimRef.cosAngle[b][2]);
  }
}

void CheckpointOutput::printRandomNumbers()
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file 
  // to read back we can use the load function
  const int N = 624;
  uint32_t* saveArray = new uint32_t[N];
  prngRef.GetGenerator()->save(saveArray);
  for(int i=0; i<N; i++) {
    outputUintIn8Chars(saveArray[i]);
  }

  // Save the location of pointer in state
  uint32_t location = prngRef.GetGenerator()->pNext -
                      prngRef.GetGenerator()->state;
  outputUintIn8Chars(location);

  // save the "left" value so we can restore it later
  outputUintIn8Chars(prngRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  outputUintIn8Chars(prngRef.GetGenerator()->seedValue);
}

void CheckpointOutput::printCoordinates()
{
  // first let's print the count
  uint32_t count = coordCurrRef.Count();
  outputUintIn8Chars(count);

  // now let's print the coordinates one by one
  for(int i=0; i<count; i++) {
    outputDoubleIn8Chars(coordCurrRef[i].x);
    outputDoubleIn8Chars(coordCurrRef[i].y);
    outputDoubleIn8Chars(coordCurrRef[i].z);
  }
}

void CheckpointOutput::printMoleculeLookupData()
{
  // print the size of molLookup array
  outputUintIn8Chars(molLookupRef.molLookupCount);
  // print the molLookup array itself
  for(int i=0; i<molLookupRef.molLookupCount; i++) {
    outputUintIn8Chars(molLookupRef.molLookup[i]);
  }

  // print the size of boxAndKindStart array
  outputUintIn8Chars(molLookupRef.boxAndKindStartCount);
  // print the BoxAndKindStart array
  for(int i=0; i<molLookupRef.boxAndKindStartCount; i++) {
    outputUintIn8Chars(molLookupRef.boxAndKindStart[i]);
  }

  // print numKinds
  outputUintIn8Chars(molLookupRef.numKinds);

  //print the size of fixedAtom array
  outputUintIn8Chars((uint)molLookupRef.fixedAtom.size());
  //print the fixedAtom array itself
  for(int i=0; i<molLookupRef.fixedAtom.size(); i++) {
    outputUintIn8Chars(molLookupRef.fixedAtom[i]);
  }
}

void CheckpointOutput::printMoveSettingsData()
{
  uint size_x, size_y, size_z;

  // print size of scale
  size_x = moveSetRef.scale.size();
  size_y = moveSetRef.scale[0].size();
  size_z = moveSetRef.scale[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print scale array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputDoubleIn8Chars(moveSetRef.scale[i][j][k]);
      }
    }
  }

  // print size of acceptPercent
  size_x = moveSetRef.acceptPercent.size();
  size_y = moveSetRef.acceptPercent[0].size();
  size_z = moveSetRef.acceptPercent[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print acceptPercent array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputDoubleIn8Chars(moveSetRef.acceptPercent[i][j][k]);
      }
    }
  }

  // print size of accepted
  size_x = moveSetRef.accepted.size();
  size_y = moveSetRef.accepted[0].size();
  size_z = moveSetRef.accepted[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print accepted array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputUintIn8Chars(moveSetRef.accepted[i][j][k]);
      }
    }
  }

  // print size of tries
  size_x = moveSetRef.tries.size();
  size_y = moveSetRef.tries[0].size();
  size_z = moveSetRef.tries[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print tries array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputUintIn8Chars(moveSetRef.tries[i][j][k]);
      }
    }
  }

  // print size of tempAccepted
  size_x = moveSetRef.tempAccepted.size();
  size_y = moveSetRef.tempAccepted[0].size();
  size_z = moveSetRef.tempAccepted[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print tempAccepted array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputUintIn8Chars(moveSetRef.tempAccepted[i][j][k]);
      }
    }
  }

  // print size of tempTries
  size_x = moveSetRef.tempTries.size();
  size_y = moveSetRef.tempTries[0].size();
  size_z = moveSetRef.tempTries[0][0].size();
  outputUintIn8Chars(size_x);
  outputUintIn8Chars(size_y);
  outputUintIn8Chars(size_z);

  // print tempTries array
  for(int i=0; i<size_x; i++) {
    for(int j=0; j<size_y; j++) {
      for(int k=0; k<size_z; k++) {
        outputUintIn8Chars(moveSetRef.tempTries[i][j][k]);
      }
    }
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

void CheckpointOutput::outputDoubleIn8Chars(double data)
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

void CheckpointOutput::outputUintIn8Chars(uint32_t data)
{
  if(outputFile == NULL) {
    fprintf(stderr, "Error opening checkpoint output file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  uint32_output_union temp;
  temp.uint_value = data;
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