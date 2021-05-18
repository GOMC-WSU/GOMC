/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <Checkpoint.h>  
  
Checkpoint::Checkpoint(){}
Checkpoint::~Checkpoint(){}
void Checkpoint::GatherCheckpointData(MoveSettings & movSetRef,
                                    MoleculeLookup & molLookRef,
                                    Molecules & molRef,
                                    PRNG & prng){
    
}
#if GOMC_LIB_MPI
void Checkpoint::GatherCheckpointData(MoveSettings & movSetRef,
                                    MoleculeLookup & molLookRef,
                                    Molecules & molRef,
                                    PRNG & prng,
                                    bool & parallelTemperingEnabled,
                                    PRNG & prngPT){
    
}
#endif
void Checkpoint::SetCheckpointData(MoveSettings & movSetRef,
                                    MoleculeLookup & molLookRef,
                                    Molecules & molRef,
                                    PRNG & prng){
    
}
#if GOMC_LIB_MPI
void Checkpoint::SetCheckpointData(MoveSettings & movSetRef,
                                    MoleculeLookup & molLookRef,
                                    Molecules & molRef,
                                    PRNG & prng,
                                    bool & parallelTemperingEnabled,
                                    PRNG & prngPT){
    
}
#endif