#ifndef CPU_SIDE_H
#define CPU_SIDE_H

//Member vars
#include "Clock.h"
#include "ConsoleOutput.h"
#include "PDBOutput.h"
#include "BlockOutput.h"
#include "FluctOutput.h"
#include "HistOutput.h"
#include "ConfigSetup.h"
#include "OutputVars.h"
#include "EnPartCntSampleOutput.h"

#include <vector>

class System;
class StaticVals;
class OutputableBase;

struct CPUSide
{

  CPUSide(System & sys, StaticVals & statV);
  void Init(PDBSetup const& pdbSet, config_setup::Output const& out,
            const ulong tillEquil, const ulong totSteps);

  void Output(const ulong step);

private:
  Clock timer;
  std::vector<OutputableBase *> outObj;
  ConsoleOutput console;
  PDBOutput pdb;
  BlockAverages block;
  Fluctuations fluct;
  Histogram hist;
#if ENSEMBLE == GCMC
  EnPartCntSample sample_N_E;
#endif

  OutputVars varRef;

};

#endif /*CPU_SIDE_H*/
