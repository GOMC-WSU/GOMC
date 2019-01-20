/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "CheckpointOutput.h"
#include "System.h"

CheckpointOutput::CheckpointOutput(System & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef),  molRef(statV.mol), 
  coordCurrRef(sys.coordinates), comCurrRef(sys.com)
{

}

void CheckpointOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output)
{
  enableOutCheckpoint = output.checkpoint.enable;
}

void CheckpointOutput::DoOutput(const ulong step)
{
  if(enableOutCheckpoint) {

  }
}