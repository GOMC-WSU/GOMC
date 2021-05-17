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

CheckpointSetup::CheckpointSetup(System & sys, StaticVals const& statV,
                                 Setup const& set) :
  parallelTemperingWasEnabled(false), molLookupRef(sys.molLookupRef), moveSetRef(sys.moveSettings)
{
  std::string file = set.config.in.files.checkpoint.name[0];
#if GOMC_LIB_MPI
  filename = sys.ms->replicaInputDirectoryPath + file;
#else
  filename = file;
#endif
  saveArray = NULL;
  #if GOMC_LIB_MPI
  saveArrayPT = NULL;
  #endif
}

CheckpointSetup::CheckpointSetup(std::string file, 
                                  MoleculeLookup & molLookup,
                                  MoveSettings & moveSettings) :
  parallelTemperingWasEnabled(false), molLookupRef(molLookup), moveSetRef(moveSettings)
{
  filename = file;
  saveArray = nullptr;
  #if GOMC_LIB_MPI
  saveArrayPT = nullptr;
  #endif
}

std::string CheckpointSetup::getFileName(){
  return filename;
}

void CheckpointSetup::SetStepNumber(ulong & startStep)
{
  startStep = stepNumber;
}

void CheckpointSetup::SetPRNGVariables(PRNG & prng)
{
  prng.GetGenerator()->load(saveArray);
  prng.GetGenerator()->pNext = prng.GetGenerator()->state + seedLocation;
  prng.GetGenerator()->left = seedLeft;
  prng.GetGenerator()->seedValue = seedValue;
}

void CheckpointSetup::SetMolecules(Molecules& mols)
{
  for(int i = 0; i < mols.count + 1; i++) {
    mols.originalStart[i] = originalStartLocalCopy[i];
  }
  for(int i = 0; i < mols.kIndexCount; i++) {
    mols.originalKIndex[i] = originalKIndexLocalCopy[i];
  }
}

bool CheckpointSetup::CheckIfParallelTemperingWasEnabled()
{
  return (bool)parallelTemperingWasEnabled;
}


#if GOMC_LIB_MPI
void CheckpointSetup::SetPRNGVariablesPT(PRNG & prng)
{
  prngPT.GetGenerator()->load(saveArrayPT);
  prngPT.GetGenerator()->pNext = prngPT.GetGenerator()->state + seedLocationPT;
  prngPT.GetGenerator()->left = seedLeftPT;
  prngPT.GetGenerator()->seedValue = seedValuePT;
}
#endif
