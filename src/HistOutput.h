/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef HIST_OUTPUT_H
#define HIST_OUTPUT_H

#include <fstream>
#include <string>

#include "../lib/BasicTypes.h" //For ulong, uint
#include "../lib/StrLib.h"
#include "EnergyTypes.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "PDBSetup.h" //For atoms class.

struct Histogram : OutputableBase {
  Histogram(OutputVars &v);

  ~Histogram();

  virtual void Sample(const ulong step);

  // No additional init.
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  void PrintKindHist(const uint b, const uint k);

  std::string GetFName(std::string const &histName, std::string const &histNum,
                       std::string const &histLetter, const uint box,
                       const uint totKinds);

  // Indices 1: boxes 2: kinds 3: count bins up to N_total
  uint **molCount[BOXES_WITH_U_NB];
  uint *total;
  uint stepsPerSample;

  std::ofstream *outF[BOXES_WITH_U_NB];
  std::string *name[BOXES_WITH_U_NB];
};

#endif /*HIST_OUTPUT_H*/
