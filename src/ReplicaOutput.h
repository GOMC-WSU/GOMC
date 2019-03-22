/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef REPLICAOUTPUT_H
#define REPLICAOUTPUT_H

#include "Writer.h"
#include "ReplicaExchange.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"


struct ReplicaOutput : OutputableBase
{
public:
  ReplicaOutput();
  ~ReplicaOutput();
  //Replica Output does not need to sample, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output)
  {
      std::string aliasStr = "";
      enableOutState = output.replicaExchange.enable;
      /*enableRestOut = output.restart.settings.enable;
      enableOut = enableOutState | enableRestOut;
      stepsCoordPerOut = output.state.settings.frequency;
      stepsRestPerOut = output.restart.settings.frequency;
      if (stepsCoordPerOut < stepsRestPerOut)
        stepsPerOut = output.state.settings.frequency;
      else
        stepsPerOut = output.restart.settings.frequency;*/
      if (enableOutState) {
          //Get alias string, based on box #.
          aliasStr = "Output Replica Log file";
          bool notify;
    #ifndef NDEBUG
          notify = true;
    #else
          notify = false;
    #endif
          outF.Init(output.state.files.replicaLog.name, aliasStr, true, notify);
          outF.open();

          recKeep.nattempt[0] = 0;
          recKeep.nattempt[1] = 0;
          recKeep.prob = (double*)malloc(output.numberOfReplicas * sizeof(double));
          recKeep.prob_sum = (double*)malloc(output.numberOfReplicas * sizeof(double));
          recKeep.nexchange = (int*)malloc(output.numberOfReplicas*sizeof(int));
          recKeep.nmoves = (int**)malloc(output.numberOfReplicas*sizeof(int*));
          for (int i = 0; i < output.numberOfReplicas; i++) {
            recKeep.nmoves[i] = (int*)malloc(output.numberOfReplicas*sizeof(int));
          }
          recKeep.indices = (int*)malloc(output.numberOfReplicas * sizeof(int));
          recKeep.p_indices = (int*)malloc(output.numberOfReplicas * sizeof(int));

          for (int i = 0; i < output.numberOfReplicas; i++){
            recKeep.indices[i] = i;
            recKeep.p_indices[i] = i;
          }
          equilSteps = output.equilSteps;
          stepsPerOut = output.stepsPerReplicaExchange;
          OutputableBase::stepsPerOut = output.stepsPerReplicaExchange;
          OutputableBase::enableOut = true;
          //DoOutput(0);
      }
  }
  virtual void DoOutput(const ulong step);
  void swapIndices(int j);
  void RecordExchangeStatistics(int j, double probability, bool bExchanged);

private:
  Writer outF;
  RecordKeeper recKeep;
  bool enableOutState;
  ulong equilSteps;
  ulong stepsPerOut;
};

#endif