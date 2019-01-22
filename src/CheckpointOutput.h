/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once

#include "OutputAbstracts.h"
#include "MoveSettings.h"
#include "Coordinates.h"
#include <iostream>

class CheckpointOutput : OutputableBase
{
public:
  CheckpointOutput(System & sys, StaticVals const& statV);

  ~CheckpointOutput();

  virtual void DoOutput(const ulong step);
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

private:
  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  Coordinates & coordCurrRef;
  COM & comCurrRef;
  PRNG & prngRef;

  bool enableOutCheckpoint;
  std::string filename;
  FILE* outputFile;

  void openOutputFile();
  void printStepNumber(const ulong step);
  void printRandomNumbers();
  void printCoordinates();
  void printMoleculeLookupData();
  void outputDoubleIn8Chars(double data);
  void outputUintIn8Chars(uint32_t data);
};
