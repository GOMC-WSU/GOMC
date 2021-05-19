/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <Checkpoint.h>  
  

Checkpoint::Checkpoint(ulong & startStep,
                        MoveSettings & movSetRef,
                        PRNG & prngRef,
                        Molecules & molRef,
                        MoleculeLookup & molLookRef){

    stepNumber = startStep;

    // Original molecule start positions.  Could be generated through kind,
    // but this allows for parallelized output.
    originalStartVec.assign(molRef.originalStart, molRef.originalStart + molRef.count + 1);
    originalKIndexVec.assign(molRef.originalKIndex, molRef.originalKIndex + molRef.kIndexCount);

    // Molecule Indices for consistent trajectories 
    originalMoleculeIndicesVec.assign(molLookRef.originalMoleculeIndices, molLookRef.originalMoleculeIndices + molLookRef.molLookupCount);
    permutedMoleculeIndicesVec.assign(molLookRef.permutedMoleculeIndices, molLookRef.permutedMoleculeIndices + molLookRef.molLookupCount);

    //ulong stepNumber;
    prngRef.GetGenerator()->save(saveArray);
    seedLocation = prngRef.GetGenerator()->pNext - prngRef.GetGenerator()->state;

    // save the "left" value so we can restore it later
    seedLeft = (prngRef.GetGenerator()->left);

    // let's save seedValue just in case
    // not sure if that is used or not, or how important it is
    seedValue = prngRef.GetGenerator()->seedValue;

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
Checkpoint::Checkpoint(){}
Checkpoint::~Checkpoint(){}
