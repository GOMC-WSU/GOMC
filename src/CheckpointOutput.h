/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once

#include "OutputAbstracts.h"
#include "MoveSettings.h"
#include "Coordinates.h"
#include "MoveBase.h"
#include <iostream>
#include "GOMC_Config.h"

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
  virtual void DoOutputRestart(const ulong step);  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}


private:
  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  PRNG & prngRef;
  Coordinates & coordCurrRef;
#if GOMC_LIB_MPI
  PRNG & prngPTRef;
#endif

  bool enableParallelTempering;
  std::string filename;
  FILE* outputFile;
  ulong stepsPerCheckpoint;
  char gomc_version[5];

  void openOutputFile();
  void setGOMCVersion();
  void printGOMCVersion();
  void printParallelTemperingBoolean();
  void printStepNumber(ulong step);
  void printRandomNumbers();
#if GOMC_LIB_MPI
  void printRandomNumbersParallelTempering();
#endif
  void printMoveSettingsData();

  void printVector3DDouble(const std::vector< std::vector< std::vector<double> > > &data);
  void printVector3DUint(const std::vector< std::vector< std::vector<uint> > > &data);
  void printVector2DUint(const std::vector< std::vector< uint > > &data);
  void printVector1DDouble(const std::vector< double > &data);
  void write_double_binary(double data);
  void write_uint8_binary(int8_t data);
  void write_uint32_binary(uint32_t data);
  void write_uint64_binary(uint64_t data);

};
