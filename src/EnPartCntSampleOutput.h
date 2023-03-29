/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef EN_PART_CNT_SAMPLE_OUTPUT_H
#define EN_PART_CNT_SAMPLE_OUTPUT_H

#include <fstream>
#include <string>

#include "BasicTypes.h" //For ulong, uint
#include "EnergyTypes.h"
#include "EnsemblePreprocessor.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "PDBSetup.h" //For atoms class.
#include "StrLib.h"

#if ENSEMBLE == GCMC

namespace config_setup {
class Output;
}

struct EnPartCntSample : OutputableBase {
  EnPartCntSample(OutputVars &v) {
    this->var = &v;
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      samplesE[b] = NULL;
      ;
      samplesN[b] = NULL;
    }
  }

  ~EnPartCntSample();

  virtual void Sample(const ulong step);

  // No additional init.
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  void WriteHeader(void);

  void InitVals(config_setup::EventSettings const &event) {
    stepsPerOut = event.frequency;
    enableOut = event.enable;
  }

  std::string GetFName(std::string const &sampleName,
                       std::string const &histNum,
                       std::string const &histLetter, const uint b);

  // samplesE --> per box; samplesN --> per kind, per box
  double *samplesE[BOXES_WITH_U_NB];
  uint **samplesN[BOXES_WITH_U_NB], stepsPerSample, samplesCollectedInFrame;
  std::ofstream outF[BOXES_WITH_U_NB];
  std::string name[BOXES_WITH_U_NB];
};

#endif /*ENSEMBLE==GCMC*/

#endif /*EN_PART_CNT_SAMPLE_OUTPUT_H*/
