/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
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

