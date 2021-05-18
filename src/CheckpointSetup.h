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

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>


class CheckpointSetup
{
public:
  CheckpointSetup();
  /* For GTesting */
  CheckpointSetup();

  ~CheckpointSetup()
  {
    if(saveArray != nullptr) {
      delete [] saveArray;
      saveArray = nullptr;
    }
    if(originalStartLocalCopy != nullptr){
      delete [] originalStartLocalCopy;
      originalStartLocalCopy = nullptr;
    }
    if(originalKIndexLocalCopy != nullptr){
      delete [] originalKIndexLocalCopy;
      originalKIndexLocalCopy = nullptr;
    }
    #if GOMC_LIB_MPI
    if(saveArrayPT != nullptr) {
      delete [] saveArrayPT;
      saveArrayPT = nullptr;
    }
    #endif
  }
  std::string getFileName();
  void SetStepNumber(ulong & startStep);
  void SetPRNGVariables(PRNG & prng);
  bool CheckIfParallelTemperingWasEnabled();
#if GOMC_LIB_MPI
  void SetPRNGVariablesPT(PRNG & prng);
#endif
  void SetMolecules(Molecules & molRef);

private:


  std::string filename;
  MoleculeLookup & molLookupRef;
  MoveSettings & moveSetRef;

  // the following variables will hold the data read from checkpoint
  // and will be passed to the rest of the code via Get functions
  bool enableParallelTemperingBool;
  int8_t parallelTemperingWasEnabled;
  char gomc_version[5];
  uint64_t stepNumber;
  // Molecules must be stored locally for InitOver
  uint32_t * originalStartLocalCopy, * originalKIndexLocalCopy;
  //ulong stepNumber;
  uint32_t* saveArray;
  uint32_t seedLocation, seedLeft, seedValue;

  // Move Settings Vectors
  std::vector<std::vector<std::vector<double> > > scaleVec, acceptPercentVec;
  std::vector<std::vector<std::vector<uint32_t> > > acceptedVec, triesVec, tempAcceptedVec,
      tempTriesVec;
  std::vector< std::vector< uint32_t > > mp_acceptedVec, mp_triesVec;
  std::vector< double > mp_r_maxVec;
  std::vector< double > mp_t_maxVec;



  #if GOMC_LIB_MPI
  uint32_t* saveArrayPT;
  uint32_t seedLocationPT, seedLeftPT, seedValuePT;
#endif

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
        assert(saveArray == nullptr);
        saveArray = new uint32_t[N + 1];
        assert(originalStartLocalCopy == nullptr);
        originalStartLocalCopy = new uint32_t[molLookupRef.molLookupCount + 1];
        assert(originalKIndexLocalCopy == nullptr);
        originalKIndexLocalCopy = new uint32_t[molLookupRef.molLookupCount];
        #if GOMC_LIB_MPI
        assert(saveArrayPT == nullptr);
        saveArrayPT = new uint32_t[N + 1];        
        #endif       
    }
    // GOMC Version
    ar & boost::serialization::make_array<char>(gomc_version, 5);  
    // Step
    ar & stepNumber;
    // PRNG Vars
    ar & boost::serialization::make_array<uint>(saveArray, N + 1);  
    ar & seedLocation;
    ar & seedLeft;
    ar & seedValue;

    // Move Settings Vectors
    ar & scaleVec;
    ar & acceptPercentVec;
    ar & acceptedVec;
    ar & triesVec;
    ar & tempAcceptedVec;
    ar & tempTriesVec;
    ar & mp_triesVec;
    ar & mp_acceptedVec;
    ar & mp_t_maxVec;
    ar & mp_r_maxVec;
    // Start and KIndex arrays
    ar & boost::serialization::make_array<uint32_t>(originalStartLocalCopy, molLookupRef.molLookupCount + 1);  
    ar & boost::serialization::make_array<uint32_t>(originalKIndexLocalCopy, molLookupRef.molLookupCount);  
    // Sorted Molecule Indices
    ar & boost::serialization::make_array<uint32_t>(molLookupRef.originalMoleculeIndices, molLookupRef.molLookupCount);  
    ar & boost::serialization::make_array<uint32_t>(molLookupRef.permutedMoleculeIndices, molLookupRef.molLookupCount);  
    // PT boolean
    ar & parallelTemperingWasEnabled;
    #if GOMC_LIB_MPI
    if((bool)parallelTemperingWasEnabled){
      // PRNG PT Vars
      ar & boost::serialization::make_array<uint>(saveArrayPT, N + 1);  
      ar & seedLocationPT;
      ar & seedLeftPT;
      ar & seedValuePT;
    }
    #endif
    std::cout << "Checkpoint loaded from " << filename << std::endl;
  }
};
