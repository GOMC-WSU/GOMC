/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <Checkpoint.h>  
  

Checkpoint::Checkpoint(const ulong & step,
                        const ulong & trueStep,
                        MoveSettings & movSetRef,
                        PRNG & prngRef,
                        const Molecules & molRef,
                        MoleculeLookup & molLookupRef,
                        MolSetup & molSetupRef,
                        pdb_setup::Atoms const& pdbSetupAtomsRef){
    GatherStep(step);
    GatherTrueStep(trueStep);
    GatherMoveSettings(movSetRef);
    GatherRandomNumbers(prngRef);
    GatherRestartMoleculeIndices(molLookupRef);
    GatherMoleculeLookup(molLookupRef);
    // Not sure if these need to be gathered..
    GatherMolSetup(molSetupRef);
    GatherPDBSetupAtoms(pdbSetupAtomsRef);
    // Not sure if these need to be gathered..
    GatherRestartMoleculeStartVec(molLookupRef, molRef);
    GatherOriginalMoleculeStartVec(molRef);
}

#if GOMC_LIB_MPI
Checkpoint::Checkpoint(const ulong & step,
                        const ulong & trueStep,
                        MoveSettings & movSetRef,
                        PRNG & prngRef,
                        const Molecules & molRef,
                        MoleculeLookup & molLookupRef,
                        bool & parallelTemperingIsEnabled,
                        PRNG & prngPTRef){
    GatherStep(step);
    GatherTrueStep(trueStep);
    GatherMoveSettings(movSetRef);
    GatherRandomNumbers(prngRef);
    GatherRestartMoleculeIndices(molLookupRef);
    GatherMoleculeLookup(molLookupRef);
    // Not sure if these need to be gathered..
    GatherMolSetup(molSetupRef);
    GatherPDBSetupAtoms(pdbSetupAtomsRef);
    // Not sure if these need to be gathered..
    GatherRestartMoleculeStartVec(molLookupRef, molRef);
    GatherOriginalMoleculeStartVec(molRef);
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


void Checkpoint::GatherStep(const ulong & step){
    stepNumber = step;
}

void Checkpoint::GatherTrueStep(const ulong & trueStep){
    trueStepNumber = trueStep;
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

/* 
   Got rid of the if condition of whether this was a restFromChk or not
   since the new method for restart from chk should believe it is always
   NOT restarting from chk. 
*/
void Checkpoint::GatherRestartMoleculeIndices(MoleculeLookup & molLookupRef){
  molLookupRef.restartMoleculeIndices.clear();
  molLookupRef.restartMoleculeIndices.resize(molLookupRef.molLookupCount);
  uint molCounter = 0, b, k, kI, countByKind, molI;
  for (b = 0; b < BOX_TOTAL; ++b) {
    for (k = 0; k < molLookupRef.numKinds; ++k) {
      countByKind = molLookupRef.NumKindInBox(k, b);
      for (kI = 0; kI < countByKind; ++kI) {
        molI = molLookupRef.GetMolNum(kI, k, b);
        molLookupRef.restartMoleculeIndices[molCounter] = molI;
        ++molCounter;
      }
    }
  }
}

void Checkpoint::GatherMoleculeLookup(MoleculeLookup & molLookupRef){
  originalMoleculeLookup = molLookupRef;
}
/* 
  After the first run, we don't parse PSF files.  We simply load the original
  molecule map from file.  Therefore, the new simulation doesn't have molecule ranges
  for the PDB data.  We generate the molecule ranges of the reordered PDB Restart
  files in this method.  
*/
void Checkpoint::GatherRestartMoleculeStartVec(MoleculeLookup & molLookupRef,
                                                const Molecules & molRef){
  restartedStartVec.clear();
  uint len = 0, sum = 0;
  //Start particle numbering @ 1
  for (uint box = 0; box < BOX_TOTAL; ++box) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(box),
                                  end = molLookupRef.BoxEnd(box);
    while (m != end) {
      uint32_t start;
      restartedStartVec.push_back(sum + len);
      molRef.GetRangeStartLength(start, len, *m);
      sum += len;
      ++m;
    }
    restartedStartVec.push_back(sum);
  }
}

void Checkpoint::GatherOriginalMoleculeStartVec(const Molecules & molRef){
  for (int i = 0; i <= molRef.count; ++i)
    originalStartVec.push_back(molRef.start[i]);
}

void Checkpoint::GatherMolSetup(MolSetup & molSetupRef){
  originalMolSetup = molSetupRef;
}

void Checkpoint::GatherPDBSetupAtoms(pdb_setup::Atoms const& pdbSetupAtomsRef){
  originalAtoms = pdbSetupAtomsRef;
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