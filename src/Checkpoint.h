/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "GOMC_Config.h"
#include <map>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/array.hpp>

#include "MoveSettings.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "PRNG.h"
#include <stdint.h>
// So we can checkpoint the MoleculeLookup
#include "MolSetup.h"
#include "PDBSetup.h"

class Checkpoint
{
    public:
    Checkpoint(const ulong & startStep,
                const ulong & trueStep,
                MoveSettings & movSetRef,
                PRNG & prng,
                const Molecules & molRef,
                MoleculeLookup & molLookRef,
                MolSetup & molSetupRef);

#if GOMC_LIB_MPI
    Checkpoint(const ulong & startStep,
                const ulong & trueStep,
                MoveSettings & movSetRef,
                PRNG & prng,
                const Molecules & molRef,
                MoleculeLookup & molLookRef,
                MolSetup & molSetupRef,
                bool & parallelTemperingIsEnabled,
                PRNG & prngPTRef);
#endif

    Checkpoint();
    ~Checkpoint();


    private:
        friend class CheckpointSetup;
        friend class CheckpointOutput;

        void GatherGOMCVersion();
        void GatherStep(const ulong & startStep);
        void GatherTrueStep(const ulong & trueStep);
        void GatherMolecules(const Molecules & molRef);
        void GatherMoleculeKindDictionary(const Molecules & molRef);
        void GatherMoveSettings(MoveSettings & movSetRef);
        void GatherSortedMoleculeIndices(MoleculeLookup & molLookupRef);
        void GatherMolSetup(MolSetup & molSetupRef);
        void GatherPDBSetupAtoms(pdb_setup::Atoms & pdbSetupAtomsRef);
        void GatherRandomNumbers(PRNG & prngRef);
    #if GOMC_LIB_MPI
        void GatherParallelTemperingBoolean(bool & parallelTemperingIsEnabled);
        void GatherRandomNumbersParallelTempering(PRNG & prngPTRef);
    #endif


        // the following variables will hold the data read from checkpoint
        // and will be passed to the rest of the code via Get functions

        char gomc_version[5];

        uint64_t stepNumber;

        // To avoid repeating Random numbers
        // on the GPU, when InitStep is set to 0
        // we maintain the true step had it not
        // been overwritten by InitStep
        // If init step isn't used
        // trueStep == step
        uint64_t trueStepNumber;

        // Original molecule start positions.  Could be generated through kind,
        // but this allows for parallelized output.
        std::vector<uint32_t> originalStartVec, originalKIndexVec;

        // Molecule Indices for consistent trajectories 
        std::vector<uint32_t> originalMoleculeIndicesVec, permutedMoleculeIndicesVec;

        // Kind indices and name map
        std::map<std::string, uint32_t> originalNameIndexMap;

        #define N_array_size 624

        //ulong stepNumber;
        uint32_t saveArray[N_array_size+1];
        uint32_t seedLocation, seedLeft, seedValue;

        // Move Settings Vectors
        std::vector<std::vector<std::vector<double> > > scaleVec, acceptPercentVec;
        std::vector<std::vector<std::vector<uint32_t> > > acceptedVec, triesVec, tempAcceptedVec,
            tempTriesVec;
        std::vector< std::vector< uint32_t > > mp_acceptedVec, mp_triesVec;
        std::vector< double > mp_r_maxVec;
        std::vector< double > mp_t_maxVec;

        #if GOMC_LIB_MPI
        int8_t parallelTemperingEnabled;

        uint32_t saveArrayPT[N_array_size+1];
        uint32_t seedLocationPT, seedLeftPT, seedValuePT;
        #endif

        MolSetup originalMolSetup;
        pdb_setup::Atoms originalAtoms;

        friend class cereal::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            // GOMC Version
            ar & gomc_version;
            // Step
            ar & stepNumber;
            ar & trueStepNumber;
            // PRNG Vars
            ar & saveArray;
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
            // name & index Map
            ar & originalNameIndexMap;
            // MolSetup
            ar & originalMolSetup;

            #if GOMC_LIB_MPI
            // PT boolean
            ar & parallelTemperingEnabled;
            if((bool)parallelTemperingEnabled){
                // PRNG PT Vars
                ar & saveArrayPT;
                ar & seedLocationPT;
                ar & seedLeftP  T;
                ar & seedValuePT;
            }
            #endif
        }
};

#endif