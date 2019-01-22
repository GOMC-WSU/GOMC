/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
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
  coordCurrRef(sys.coordinates), comCurrRef(sys.com),
  filename("checkpoint.dat")
{

}

void CheckpointOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output)
{
  enableOutCheckpoint = output.checkpoint.enable;
}

void CheckpointOutput::DoOutput(const ulong step)
{
  if(enableOutCheckpoint) {
    openOutputFile();
    printStepNumber(step);
    printRandomNumbers();
    printCoordinates();
    printMoleculeLookupData();
  }
}

void CheckpointOutput::printStepNumber(const ulong step)
{
  uint32_t s = step;
  outputUintIn8Chars(s);
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

  // print the canSwapKind array size
  outputUintIn8Chars(molLookupRef.canSwapKind.size());

  //print the canSwapKind array
  for(int i=0; i<molLookupRef.canSwapKind.size(); i++) {
    outputUintIn8Chars(molLookupRef.canSwapKind[i]);
  }

  // print the size of canMoveKind array
  outputUintIn8Chars(molLookupRef.canMoveKind.size());

  // print the canMoveKind array
  for(int i=0; i<molLookupRef.canMoveKind.size(); i++) {
    outputUintIn8Chars(molLookupRef.canMoveKind[i]);
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
}

void CheckpointOutput::outputUintIn8Chars(uint32_t data)
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
}