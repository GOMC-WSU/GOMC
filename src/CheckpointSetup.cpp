/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "CheckpointSetup.h"

CheckpointSetup::CheckpointSetup(ulong & startStep,
                                ulong & trueStep,
                                MoleculeLookup & molLookup, 
                                MoveSettings & moveSettings,
                                Molecules & mol,
                                PRNG & prng,
                                Random123Wrapper & r123,
                                Setup const& set) :
  molLookupRef(molLookup), moveSetRef(moveSettings), molRef(mol), prngRef(prng),
  r123Ref(r123), startStepRef(startStep), trueStepRef(trueStep)
{
  std::string file = set.config.in.files.checkpoint.name[0];
#if GOMC_LIB_MPI
  filename = sys.ms->replicaInputDirectoryPath + file;
#else
  filename = file;
#endif
}

CheckpointSetup::~CheckpointSetup(){}


std::string CheckpointSetup::getFileName(){
  return filename;
}

void CheckpointSetup::loadCheckpointFile(){
  // create and open a character archive for intput
  std::ifstream ifs(filename);
  if (!ifs.is_open()){
    fprintf(stderr, "Error opening checkpoint input file %s\n",
            filename.c_str());
    exit(EXIT_FAILURE);
  }
  cereal::BinaryInputArchive ia(ifs);
  ia >> chkObj;
  SetCheckpointData();
  std::cout << "Checkpoint loaded from " << filename << std::endl;
}

void CheckpointSetup::InitOver(Molecules & molRef){
  SetMolecules(molRef);
  SetMoleculeKindDictionary(molRef);
}

void CheckpointSetup::SetCheckpointData   (){
  SetStepNumber();
  SetTrueStepNumber();
  SetMoveSettings();
  SetPRNGVariables();
  SetR123Variables();
  SetMolecules();
  SetMoleculeKindDictionary();
  SetMoleculeIndices();
}

#if GOMC_LIB_MPI
void CheckpointSetup::SetCheckpointData   (bool & parallelTemperingIsEnabled,
                                          PRNG & prngPT){
  SetStepNumber();
  SetTrueStepNumber();
  SetMoveSettings();
  SetPRNGVariables();
  SetR123Variables();
  SetMolecules();
  SetMoleculeKindDictionary();
  SetMoleculeIndices();
  SetParallelTemperingWasEnabled();
  if(parallelTemperingIsEnabled && parallelTemperingWasEnabled)
    SetPRNGVariablesPT();
}
#endif

void CheckpointSetup::SetStepNumber()
{
  startStepRef = chkObj.stepNumber;
}

void CheckpointSetup::SetTrueStepNumber()
{
  printf("%-40s %-lu \n", "Info: Loading true step from checkpoint", chkObj.trueStepNumber);
  trueStepRef = chkObj.trueStepNumber;
}

void CheckpointSetup::SetMolecules(Molecules& mols)
{
  /* Original Start Indices are for space demarcation in trajectory frame */
  mols.originalStart = vect::transfer<uint32_t>(chkObj.originalStartVec);
  /* Kinds store accessory molecule data such as residue, charge, etc */
  mols.originalKIndex = vect::transfer<uint32_t>(chkObj.originalKIndexVec);
}

// Pass mol reference because the old reference might break after
// the InitOver call..
void CheckpointSetup::SetMoleculeKindDictionary(Molecules& mols){
  // Kind indices and name map
  std::map<std::string, uint32_t> nameIndexMap;
  for (uint mk = 0 ; mk < mols.kindsCount; mk++)
    nameIndexMap[mols.kinds[mk].name] = (uint32_t)mk;

  for (uint mk = 0 ; mk < mols.kindsCount; mk++)
    mols.originalKIndex2CurrentKIndex[chkObj.originalNameIndexMap[mols.kinds[mk].name]]
         = nameIndexMap[mols.kinds[mk].name];
}

void CheckpointSetup::SetMoleculeKindDictionary(){
  // Kind indices and name map
  std::map<std::string, uint32_t> nameIndexMap;
  for (uint mk = 0 ; mk < molRef.kindsCount; mk++)
    nameIndexMap[molRef.kinds[mk].name] = (uint32_t)mk;

  for (uint mk = 0 ; mk < molRef.kindsCount; mk++)
    molRef.originalKIndex2CurrentKIndex[chkObj.originalNameIndexMap[molRef.kinds[mk].name]]
         = nameIndexMap[molRef.kinds[mk].name];
}

void CheckpointSetup::SetMoveSettings()
{
  moveSetRef.scale = chkObj.scaleVec;
  moveSetRef.acceptPercent = chkObj.acceptPercentVec;
  moveSetRef.accepted = chkObj.acceptedVec;
  moveSetRef.tries = chkObj.triesVec;
  moveSetRef.tempAccepted = chkObj.tempAcceptedVec;
  moveSetRef.tempTries = chkObj.tempTriesVec;
  moveSetRef.mp_tries = chkObj.mp_triesVec;
  moveSetRef.mp_accepted = chkObj.mp_acceptedVec;
  moveSetRef.mp_t_max = chkObj.mp_t_maxVec;
  moveSetRef.mp_r_max = chkObj.mp_r_maxVec;
}

void CheckpointSetup::SetPRNGVariables()
{
  prngRef.GetGenerator()->load(chkObj.saveArray);
  prngRef.GetGenerator()->pNext = prngRef.GetGenerator()->state + chkObj.seedLocation;
  prngRef.GetGenerator()->left = chkObj.seedLeft;
  prngRef.GetGenerator()->seedValue = chkObj.seedValue;
}

void CheckpointSetup::SetR123Variables()
{
  r123Ref.SetRandomSeed(chkObj.seedValue);
}

void CheckpointSetup::SetMolecules()
{
  /* Original Start Indices are for space demarcation in trajectory frame */
  molRef.originalStart = vect::transfer<uint32_t>(chkObj.originalStartVec);
  /* Kinds store accessory molecule data such as residue, charge, etc */
  molRef.originalKIndex = vect::transfer<uint32_t>(chkObj.originalKIndexVec);
}

void CheckpointSetup::SetMoleculeIndices(){
  /* Original Mol Indices are for constant trajectory output from start to finish of a single run*/
  molLookupRef.originalMoleculeIndices = vect::transfer<uint32_t>(chkObj.originalMoleculeIndicesVec);
  /* Permuted Mol Indices are for following single molecules as molLookup permutes the indices and continuing the next run*/
  molLookupRef.permutedMoleculeIndices = vect::transfer<uint32_t>(chkObj.permutedMoleculeIndicesVec);
}


#if GOMC_LIB_MPI
bool CheckpointSetup::SetParallelTemperingWasEnabled()
{
  parallelTemperingWasEnabled = (bool)chkObj.parallelTemperingEnabled;
}

void CheckpointSetup::SetPRNGVariablesPT(PRNG & prng)
{
  prngPT.GetGenerator()->load(chkObj.saveArrayPT);
  prngPT.GetGenerator()->pNext = prngPT.GetGenerator()->state + chkObj.seedLocationPT;
  prngPT.GetGenerator()->left = chkObj.seedLeftPT;
  prngPT.GetGenerator()->seedValue = chkObj.seedValuePT;
}
#endif
