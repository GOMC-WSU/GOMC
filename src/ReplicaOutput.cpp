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