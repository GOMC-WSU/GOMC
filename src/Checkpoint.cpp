/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/

#include <Checkpoint.h>

Checkpoint::Checkpoint(const ulong &step, const ulong &trueStep,
                       MoveSettings &movSetRef, PRNG &prngRef,
                       const Molecules &molRef, MoleculeLookup &molLookupRef,
                       MolSetup &molSetupRef,
                       pdb_setup::Atoms const &pdbSetupAtomsRef) {
  GatherStep(step);
  GatherTrueStep(trueStep);
  GatherMoveSettings(movSetRef);
  GatherRandomNumbers(prngRef);
  GatherMoleculeLookup(molLookupRef, molRef);
  // Not sure if these need to be gathered..
  GatherMolSetup(molSetupRef);
  GatherPDBSetupAtoms(pdbSetupAtomsRef);
  // Not sure if these need to be gathered..
  GatherRestartMoleculeStartVec(molLookupRef, molRef);
  GatherOriginalMoleculeStartVec(molRef);
}

#if GOMC_LIB_MPI
Checkpoint::Checkpoint(const ulong &step, const ulong &trueStep,
                       MoveSettings &movSetRef, PRNG &prngRef,
                       const Molecules &molRef, MoleculeLookup &molLookupRef,
                       bool &parallelTemperingIsEnabled, PRNG &prngPTRef) {
  GatherStep(step);
  GatherTrueStep(trueStep);
  GatherMoveSettings(movSetRef);
  GatherRandomNumbers(prngRef);
  GatherMoleculeLookup(molLookupRef, molRef);
  // Not sure if these need to be gathered..
  GatherMolSetup(molSetupRef);
  GatherPDBSetupAtoms(pdbSetupAtomsRef);
  // Not sure if these need to be gathered..
  GatherRestartMoleculeStartVec(molLookupRef, molRef);
  GatherOriginalMoleculeStartVec(molRef);
  GatherParallelTemperingBoolean(parallelTemperingIsEnabled);
  if (parallelTemperingIsEnabled)
    GatherRandomNumbersParallelTempering(prngPTRef);
}
#endif

Checkpoint::Checkpoint() {}
Checkpoint::~Checkpoint() {}

void Checkpoint::GatherGOMCVersion() {
  sprintf(gomc_version, "%d.%02d", GOMC_VERSION_MAJOR,
          GOMC_VERSION_MINOR % 100);
}

void Checkpoint::GatherStep(const ulong &step) { stepNumber = step; }

void Checkpoint::GatherTrueStep(const ulong &trueStep) {
  trueStepNumber = trueStep;
}

void Checkpoint::GatherRandomNumbers(PRNG &prngRef) {
  prngRef.GetGenerator()->save(saveArray);
  seedLocation = prngRef.GetGenerator()->pNext - prngRef.GetGenerator()->state;

  // save the "left" value so we can restore it later
  seedLeft = (prngRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedValue = (prngRef.GetGenerator()->seedValue);
}

/* After the first run, the molecules are sorted, so we need to use the same
   sorting process seen below, to reinitialize the originalMolInds every
   checkpoint */
void Checkpoint::GatherMoveSettings(MoveSettings &movSetRef) {
  // Move Settings Vectors
  scaleVec = movSetRef.scale;
  acceptPercentVec = movSetRef.acceptPercent;
  acceptedVec = movSetRef.accepted;
  triesVec = movSetRef.tries;
  tempAcceptedVec = movSetRef.tempAccepted;
  tempTriesVec = movSetRef.tempTries;
  mp_triesVec = movSetRef.mp_tries;
  mp_acceptedVec = movSetRef.mp_accepted;
  mp_interval_acceptedVec = movSetRef.mp_interval_accepted;
  mp_interval_triesVec = movSetRef.mp_interval_tries;
  mp_r_maxVec = movSetRef.mp_r_max;
  mp_t_maxVec = movSetRef.mp_t_max;
  isSingleMoveAcceptedVec = movSetRef.isSingleMoveAccepted;
}

/* Create vector versions of the arrays for simple serialization.
   Also create a permuted list of indices indicating the load position
   of the restarted pdb/psf molecules into the original PDBAtoms data structure.
 */
void Checkpoint::GatherMoleculeLookup(MoleculeLookup &molLookupRef,
                                      const Molecules &molRef) {
  originalMoleculeLookup = molLookupRef;
}
/*
  After the first run, we don't parse PSF files.  We simply load the original
  molecule map from file.  Therefore, the new simulation doesn't have molecule
  ranges for the PDB data.  We generate the molecule ranges of the reordered PDB
  Restart files in this method.
*/
void Checkpoint::GatherRestartMoleculeStartVec(MoleculeLookup &molLookupRef,
                                               const Molecules &molRef) {
  restartedStartVec.clear();
  uint len, sum = 0, start;
  // Start particle numbering @ 1
  restartedStartVec.push_back(0);
  for (uint box = 0; box < BOX_TOTAL; ++box) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(box),
                                 end = molLookupRef.BoxEnd(box);
    while (m != end) {
      start = 0;
      len = 0;
      molRef.GetRangeStartLength(start, len, *m);
      sum += len;
      restartedStartVec.push_back(sum);
      ++m;
    }
  }
}

void Checkpoint::GatherOriginalMoleculeStartVec(const Molecules &molRef) {
  for (int i = 0; i <= molRef.count; ++i)
    originalStartVec.push_back(molRef.start[i]);
}

void Checkpoint::GatherMolSetup(MolSetup &molSetupRef) {
  originalMolSetup = molSetupRef;
}

void Checkpoint::GatherPDBSetupAtoms(pdb_setup::Atoms const &pdbSetupAtomsRef) {
  originalAtoms = pdbSetupAtomsRef;
}

#if GOMC_LIB_MPI
void Checkpoint::GatherParallelTemperingBoolean(
    bool &parallelTemperingIsEnabled) {
  parallelTemperingEnabled = parallelTemperingIsEnabled;
}

void Checkpoint::GatherRandomNumbersParallelTempering(PRNG &prngPTRef) {
  // First let's save the state array inside prng
  // the length of the array is 624
  // there is a save function inside MersenneTwister.h file
  // to read back we can use the load function
  // The MT::save function also appends the "left" variable,
  // so need to allocate one more array element
  prngPTRef.GetGenerator()->save(saveArrayPT);
  // Save the location of pointer in state
  locationPT =
      prngPTRef.GetGenerator()->pNext - prngPTRef.GetGenerator()->state;

  // save the "left" value so we can restore it later
  leftPT = (prngPTRef.GetGenerator()->left);

  // let's save seedValue just in case
  // not sure if that is used or not, or how important it is
  seedPT = prngPTRef.GetGenerator()->seedValue;
}
#endif