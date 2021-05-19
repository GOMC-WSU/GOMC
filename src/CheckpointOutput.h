/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef CHECKPOINT_OUTPUT_H
#define CHECKPOINT_OUTPUT__H

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
  }

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}


private:
  void GatherCheckpointData(MoveSettings & movSetRef,
                            MoleculeLookup & molLookRef,
                            Molecules & molRef,
                            PRNG & prng);
  #if GOMC_LIB_MPI
  void GatherCheckpointData(MoveSettings & movSetRef,
                            MoleculeLookup & molLookRef,
                            Molecules & molRef,
                            PRNG & prng,
                            bool & parallelTemperingEnabled,
                            PRNG & prngPT);
  #endif


  std::ofstream * ofs;
  boost::archive::text_oarchive * oa;

  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  PRNG & prngRef;
  Coordinates & coordCurrRef;
#if GOMC_LIB_MPI
  PRNG & prngPTRef;
#endif

  bool enableParallelTemperingBool;
  std::string filename;
  ulong stepsPerCheckpoint;

  void openOutputFile(std::string filenameArg);
  void closeOutputFile();
  void setGOMCVersion();
  void setParallelTemperingBoolean();
  void setStepNumber(ulong stepArg);
  void setRandomNumbers();
  void setSortedMoleculeIndices();

#if GOMC_LIB_MPI
  void setRandomNumbersParallelTempering();
#endif
  void setMoveSettingsData();


};

#endif
