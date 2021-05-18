/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once

#include "GOMC_Config.h"

#include "MoveSettings.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "PRNG.h"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>

class Checkpoint
{
    public:

        Checkpoint();
        ~Checkpoint();
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
        void SetCheckpointData   (MoveSettings & movSetRef,
                                  MoleculeLookup & molLookRef,
                                  Molecules & molRef,
                                  PRNG & prng);
        #if GOMC_LIB_MPI
        void SetCheckpointData   (MoveSettings & movSetRef,
                                  MoleculeLookup & molLookRef,
                                  Molecules & molRef,
                                  PRNG & prng,
                                  bool & parallelTemperingEnabled,
                                  PRNG & prngPT);
        #endif

    private:

        // the following variables will hold the data read from checkpoint
        // and will be passed to the rest of the code via Get functions

        char gomc_version[5];
        uint64_t stepNumber;

        // Original molecule start positions.  Could be generated through kind,
        // but this allows for parallelized output.
        std::vector<uint32_t> originalStartVec, originalKIndexVec;

        // Molecule Indices for consistent trajectories 
        std::vector<uint32_t> originalMoleculeIndicesVec, permutedMoleculeIndicesVec;

        #define N 624

        //ulong stepNumber;
        uint32_t saveArray[N+1];
        uint32_t seedLocation, seedLeft, seedValue;

        // Move Settings Vectors
        std::vector<std::vector<std::vector<double> > > scaleVec, acceptPercentVec;
        std::vector<std::vector<std::vector<uint32_t> > > acceptedVec, triesVec, tempAcceptedVec,
            tempTriesVec;
        std::vector< std::vector< uint32_t > > mp_acceptedVec, mp_triesVec;
        std::vector< double > mp_r_maxVec;
        std::vector< double > mp_t_maxVec;

        int8_t parallelTemperingEnabled;

        #if GOMC_LIB_MPI
        uint32_t saveArrayPT[N+1];
        uint32_t seedLocationPT, seedLeftPT, seedValuePT;
        #endif

        friend class boost::serialization::access;
        // When the class Archive corresponds to an output archive, the
        // & operator is defined similar to <<.  Likewise, when the class Archive
        // is a type of input archive the & operator is defined similar to >>.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
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
            ar & originalStartVec;  
            ar & originalKIndexVec;  
            // Sorted Molecule Indices
            ar & originalMoleculeIndicesVec;  
            ar & permutedMoleculeIndicesVec;  
            // PT boolean
            ar & parallelTemperingEnabled;
            #if GOMC_LIB_MPI
            if((bool)parallelTemperingEnabled){
                // PRNG PT Vars
                ar & boost::serialization::make_array<uint>(saveArrayPT, N + 1);  
                ar & seedLocationPT;
                ar & seedLeftPT;
                ar & seedValuePT;
            }
            #endif
        }
}