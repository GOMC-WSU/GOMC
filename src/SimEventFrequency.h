/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef SIM_EVENT_FREQUENCY_H
#define SIM_EVENT_FREQUENCY_H

#include "BasicTypes.h"  //for ulong
#include "ConfigSetup.h" //for event frequencies from config file.

struct SimEventFrequency {
  ulong total, perAdjust, tillEquil, pCalcFreq, parallelTempFreq;
  bool pressureCalc, parallelTemp;

  void Init(config_setup::Step const &s) {
    total = s.total;
    perAdjust = s.adjustment;
    tillEquil = s.equil;
    pCalcFreq = s.pressureCalcFreq;
    pressureCalc = s.pressureCalc;
#if GOMC_LIB_MPI
    parallelTempFreq = s.parallelTempFreq;
    parallelTemp = s.parallelTemp;
#endif
  }
};

#endif /*SIM_EVENT_FREQUENCY_H*/
