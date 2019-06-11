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

class CheckpointOutput : public OutputableBase
{
public:
  CheckpointOutput(System & sys, StaticVals const& statV);

  ~CheckpointOutput()
  {
    if(outputFile)
      fclose(outputFile);
  }

  virtual void DoOutput(const ulong step);
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}
  virtual void Output(const ulong step)
  {
    if(!enableOutCheckpoint) {
      return;
    }

    if((step + 1) % stepsPerCheckpoint == 0) {
      DoOutput(step);
    }
  }

private:
  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  PRNG & prngRef;
  Coordinates & coordCurrRef;

  bool enableOutCheckpoint;
  std::string filename;
  FILE* outputFile;
  ulong stepsPerCheckpoint;

  void openOutputFile();
  void printStepNumber(const ulong step);
  void printRandomNumbers();
  void printCoordinates();
  void printMoleculeLookupData();
  void printMoveSettingsData();
  void printBoxDimensionsData();

  void outputRealIn8Chars(real data);
  void outputUintIn8Chars(uint32_t data);
};
