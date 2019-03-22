/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ReplicaOutput.h"

ReplicaOutput::ReplicaOutput(){}
ReplicaOutput::~ReplicaOutput(){
    outF.close();
}
//void ReplicaOutput::Init(){}

void ReplicaOutput::DoOutput(const ulong step){
    if (step >= equilSteps){
     outF.file << "Replica exchange at step " << step << std::endl;
    }
}

void ReplicaOutput::RecordExchangeStatistics(int j, double probability, bool bExchanged){
    recKeep.prob_sum[recKeep.indices[j]] += probability;
    if (bExchanged)
        recKeep.nexchange[recKeep.indices[j]]++;
    
}

void ReplicaOutput::swapIndices(int j){
    int temp;
    temp = recKeep.indices[j+1];
    recKeep.indices[j+1] = recKeep.indices[j];
    recKeep.indices[j] = temp;
}