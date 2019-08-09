/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once

#include "OutputAbstracts.h"
#include "MoveSettings.h"
#include "Coordinates.h"
#include <iostream>

class CheckpointSetup
{
public:
  CheckpointSetup(System & sys, StaticVals const& statV);

  ~CheckpointSetup()
  {
    if(inputFile != NULL) {
      fclose(inputFile);
      inputFile = NULL;
    }
    if(saveArray != NULL) {
      delete [] saveArray;
      saveArray = NULL;
    }
  }

  void ReadAll();
  void SetStepNumber(ulong & startStep);
  void SetPRNGVariables(PRNG & prng);
  void SetBoxDimensions(BoxDimensions & boxDimRef);
  void SetCoordinates(Coordinates & coordinates);
  void SetMoleculeLookup(MoleculeLookup & molLookupRef);
  void SetMoveSettings(MoveSettings & moveSettings);

private:
  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  Coordinates & coordCurrRef;
  PRNG & prngRef;

  std::string filename;
  FILE* inputFile;

  // the following variables will hold the data read from checkpoint
  // and will be passed to the rest of the code via Get functions
  ulong stepNumber;
  uint32_t totalBoxes;
  vector<vector<double> > axis;
  vector<vector<double> > cosAngle;
  uint32_t* saveArray;
  uint32_t seedLocation, seedLeft, seedValue;
  uint32_t coordLength;
  XYZArray coords;
  vector<uint32_t> molLookupVec, boxAndKindStartVec, fixedAtomVec;
  uint32_t numKinds;
  vector<vector<vector<double> > > scaleVec, acceptPercentVec;
  vector<vector<vector<uint32_t> > > acceptedVec, triesVec, tempAcceptedVec,
         tempTriesVec;

  // private functions used by ReadAll and Get functions
  void openInputFile();
  void readStepNumber();
  void readRandomNumbers();
  void readCoordinates();
  void readMoleculeLookupData();
  void readMoveSettingsData();
  void readBoxDimensionsData();
  void closeInputFile();

  double readDoubleIn8Chars();
  uint32_t readUintIn8Chars();
};
