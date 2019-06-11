/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EN_PART_CNT_SAMPLE_OUTPUT_H
#define EN_PART_CNT_SAMPLE_OUTPUT_H

#include "EnsemblePreprocessor.h"

#include <string>
#include <fstream>

#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "BasicTypes.h" //For ulong, uint
#include "StrLib.h"
#include "PDBSetup.h" //For atoms class.
#include "EnergyTypes.h"

#if ENSEMBLE == GCMC

namespace config_setup
{
class Output;
}

struct EnPartCntSample : OutputableBase {
  EnPartCntSample(OutputVars & v)
  {
    this->var = &v;
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      samplesE[b] = NULL;;
      samplesN[b] = NULL;
    }

  }

  ~EnPartCntSample();

  virtual void Sample(const ulong step);

  //No additional init.
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void DoOutput(const ulong step);

private:
  void WriteHeader(void);

  void InitVals(config_setup::EventSettings const& event)
  {
    stepsPerOut = event.frequency;
    enableOut = event.enable;
  }

  std::string GetFName(std::string const& sampleName,
                       std::string const& histNum,
                       std::string const& histLetter,
                       const uint b);

  //samplesE --> per box; samplesN --> per kind, per box
  real * samplesE [BOXES_WITH_U_NB];
  uint ** samplesN [BOXES_WITH_U_NB], stepsPerSample, samplesCollectedInFrame;
  std::ofstream outF[BOXES_WITH_U_NB];
  std::string name [BOXES_WITH_U_NB];
};

#endif /*ENSEMBLE==GCMC*/

#endif /*EN_PART_CNT_SAMPLE_OUTPUT_H*/
