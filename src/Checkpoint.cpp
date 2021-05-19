/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <Checkpoint.h>  
  

Checkpoint::Checkpoint(const ulong & startStep,
                        MoveSettings & movSetRef,
                        PRNG & prngRef,
                        const Molecules & molRef,
                        MoleculeLookup & molLookupRef){
    GatherStep(startStep);
    GatherMoveSettings(movSetRef);
    GatherRandomNumbers(prngRef);
    GatherMolecules(molRef);
    GatherSortedMoleculeIndices(molLookupRef);
}

#if GOMC_LIB_MPI
Checkpoint::Checkpoint(const ulong & startStep,
                        MoveSettings & movSetRef,
                        PRNG & prngRef,
                        const Molecules & molRef,
                        MoleculeLookup & molLookupRef,
                        bool & parallelTemperingIsEnabled,
                        PRNG & prngPTRef){
    GatherStep(startStep);
    GatherMoveSettings(movSetRef);
    GatherRandomNumbers(prngRef);
    GatherMolecules(molRef);
    GatherSortedMoleculeIndices(molLookupRef);
    GatherParallelTemperingBoolean(parallelTemperingIsEnabled);
    if(parallelTemperingIsEnabled)
        GatherRandomNumbersParallelTempering(prngPTRef);
}
#endif

Checkpoint::Checkpoint(){}
Checkpoint::~Checkpoint(){}

void Checkpoint::GatherGOMCVersion()
{
  sprintf(gomc_version, "%d.%02d", GOMC_VERSION_MAJOR, GOMC_VERSION_MINOR % 100);
}


void Checkpoint::GatherStep(const ulong & startStep){
    stepNumber = startStep;
}


void Checkpoint::GatherMolecules(const Molecules & molRef){
    // Original molecule start positions.  Could be generated through kind,
    // but this allows for parallelized output.
    originalStartVec.assign(molRef.originalStart, molRef.originalStart + molRef.count + 1);
    originalKIndexVec.assign(molRef.originalKIndex, molRef.originalKIndex + molRef.kIndexCount);
}

void Checkpoint::GatherRandomNumbers(PRNG & prngRef){

    prngRef.GetGenerator()->save(saveArray);
    seedLocation = prngRef.GetGenerator()->pNext - prngRef.GetGenerator()->state;

    // save the "left" value so we can restore it later
    seedLeft = (prngRef.GetGenerator()->left);

    // let's save seedValue just in case
    // not sure if that is used or not, or how important it is
    seedValue = prngRef.GetGenerator()->seedValue;
}

/* After the first run, the molecules are sorted, so we need to use the same sorting process
   seen below, to reinitialize the originalMolInds every checkpoint */
void Checkpoint::GatherMoveSettings(MoveSettings & movSetRef){
    // Move Settings Vectors
    scaleVec = movSetRef.scale;
    acceptPercentVec = movSetRef.acceptPercent;
    
    acceptedVec = movSetRef.accepted;
    triesVec = movSetRef.tries;
    tempAcceptedVec = movSetRef.tempAccepted;
    tempTriesVec = movSetRef.tempTries;
    mp_acceptedVec = movSetRef.mp_accepted;
    mp_triesVec = movSetRef.mp_tries;
    mp_r_maxVec = movSetRef.mp_r_max;
    mp_t_maxVec = movSetRef.mp_t_max;
}

/* After the first run, the molecules are sorted, so we need to use the same sorting process
   seen below, to reinitialize the originalMolInds every checkpoint */
void Checkpoint::GatherSortedMoleculeIndices(MoleculeLookup & molLookupRef){
  originalMoleculeIndicesVec.clear();
  permutedMoleculeIndicesVec.clear();
  originalMoleculeIndicesVec.resize(molLookupRef.molLookupCount);
  permutedMoleculeIndicesVec.resize(molLookupRef.molLookupCount);
  uint molCounter = 0, b, k, kI, countByKind, molI;
  if (!molLookupRef.restartFromCheckpoint){
    for (b = 0; b < BOX_TOTAL; ++b) {
      for (k = 0; k < molLookupRef.numKinds; ++k) {
        countByKind = molLookupRef.NumKindInBox(k, b);
        for (kI = 0; kI < countByKind; ++kI) {
          molI = molLookupRef.GetMolNum(kI, k, b);
          originalMoleculeIndicesVec[molCounter] = molI;
          ++molCounter;
        }
      }
    }
    for (uint molI = 0; molI < molLookupRef.molLookupCount; ++molI){
      permutedMoleculeIndicesVec[molI] = molLookupRef.originalMoleculeIndices[molI];
    }
  } else {
    for (b = 0; b < BOX_TOTAL; ++b) {
      for (k = 0; k < molLookupRef.numKinds; ++k) {
        countByKind = molLookupRef.NumKindInBox(k, b);
        for (kI = 0; kI < countByKind; ++kI) {
          molI = molLookupRef.GetSortedMolNum(kI, k, b);
          originalMoleculeIndicesVec[molCounter] = molLookupRef.permutedMoleculeIndices[molI];
          ++molCounter;
        }
      }
    }
  }
}

#if GOMC_LIB_MPI
void Checkpoint::GatherParallelTemperingBoolean(bool & parallelTemperingIsEnabled){
    parallelTemperingEnabled = parallelTemperingIsEnabled;
}


void Checkpoint::GatherRandomNumbersParallelTempering(PRNG & prngPTRef)
{
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  // The MT::save function also appends the "left" variable,
  // so need to allocate one more array element
  prngPTRef.GetGenerator()->save(saveArrayPT);
  // Save the location of pointer in state
  locationPT = prngPTRef.GetGenerator()->pNext - prngPTRef.GetGenerator()->state;

  // save the "left" value so we can restore it later
  leftPT = (prngPTRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedPT = prngPTRef.GetGenerator()->seedValue;
}
#endif