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
#include "PDBSetup.h"
#include <iostream>
#include "VectorLib.h" //for transfer.


class CheckpointSetup
{
public:
  CheckpointSetup(System & sys, StaticVals const& statV, Setup const& set);
  /* For GTesting */
  CheckpointSetup(std::string file);

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
  bool CheckIfParallelTemperingWasEnabled();
#if GOMC_LIB_MPI
  void SetPRNGVariablesPT(PRNG & prng);
#endif
  void SetMoveSettings(MoveSettings & moveSettings);
  void SetMoleculeLookup(MoleculeLookup & molLookupRef);
  void SetMolecules(Molecules & molRef);
  void SetOriginalMoleculeIndices(MoleculeLookup & molLookupRef);


  // molecules data, for GTest I need this public
  // Will add GTest macros around this when I merge PTTesters Branch
  std::vector< uint > originalMoleculeIndicesVec;
  std::vector< uint > permutedMoleculeIndicesVec; 
  std::vector< uint > molecules_originalStartVec, molecules_originalKIndexVec;

private:
  std::string filename;
  FILE* inputFile;

  // the following variables will hold the data read from checkpoint
  // and will be passed to the rest of the code via Get functions
  int8_t parallelTemperingWasEnabled;
  char gomc_version[5];
  ulong stepNumber;
  uint32_t totalBoxes;
  uint32_t* saveArray;
  uint32_t seedLocation, seedLeft, seedValue;
  std::vector<uint32_t> molLookupVec, boxAndKindStartVec, fixedMoleculeVec;
  uint32_t numKinds;
  std::vector<std::vector<std::vector<double> > > scaleVec, acceptPercentVec;
  std::vector<std::vector<std::vector<uint32_t> > > acceptedVec, triesVec, tempAcceptedVec,
      tempTriesVec;
  std::vector< std::vector< uint > > mp_acceptedVec, mp_triesVec;
  std::vector< double > mp_r_maxVec;
  std::vector< double > mp_t_maxVec;



  // private functions used by ReadAll and Get functions
  void readGOMCVersion();
  bool isLegacy();
  void openInputFile();
  void readParallelTemperingBoolean();
  void readStepNumber();
  void readOriginalMoleculeIndices();
  void readRandomNumbers();
#if GOMC_LIB_MPI
  void readRandomNumbersParallelTempering();
  uint32_t* saveArrayPT;
  uint32_t seedLocationPT, seedLeftPT, seedValuePT;
#endif
  void readMoleculesData();
  void readMoleculeLookupData();
  void readMoveSettingsData();
  void closeInputFile();

  void readVector3DDouble(std::vector< std::vector< std::vector <double> > > & data);
  void readVector3DUint(std::vector< std::vector< std::vector <uint> > > & data);
  void readVector2DUint(std::vector< std::vector< uint > > & data);
  void readVector1DUint(std::vector< uint > & data);
  void readVector1DDouble(std::vector< double > & data);
  double read_double_binary();
  int8_t read_uint8_binary();
  uint32_t read_uint32_binary();
  uint64_t read_uint64_binary();

  const int N = 624;

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    if (Archive::is_loading::value)
    {
        assert(gomc_version == nullptr);
        gomc_version = new char[5];
        assert(atomKind == nullptr);
        saveArray = new uint32_t[N + 1];
        assert(atomCharge == nullptr);
        atomCharge = new double[numAtoms];        
    }
    // GOMC Version
    ar & boost::serialization::make_array<char>(gomc_version, 5);  
    // Step
    ar & step;
    // PRNG Vars
    ar & boost::serialization::make_array<uint>(saveArray, N + 1);  
    ar & location;
    ar & left;
    ar & seed;
    // Move Settings Vectors
    ar & moveSetRef.scale;
    ar & moveSetRef.acceptPercent;
    ar & moveSetRef.accepted;
    ar & moveSetRef.tries;
    ar & moveSetRef.tempAccepted;
    ar & moveSetRef.tempTries;
    ar & moveSetRef.mp_tries;
    ar & moveSetRef.mp_accepted;
    ar & moveSetRef.mp_t_max;
    ar & moveSetRef.mp_r_max;
    // Sorted Molecule Indices
    ar & originalMoleculeIndicesVec;
    ar & boost::serialization::make_array<uint>(molLookupRef.permutedMoleculeIndices, molLookupRef.molLookupCount);  
    // PT boolean
    ar & enableParallelTempering;
    #if GOMC_LIB_MPI
      // PRNG PT Vars
      ar & boost::serialization::make_array<uint>(saveArrayPT, N + 1);  
      ar & locationPT;
      ar & leftPT;
      ar & seedPT;
    #endif
  }
};
