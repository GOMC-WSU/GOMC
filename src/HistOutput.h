/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef HIST_OUTPUT_H
#define HIST_OUTPUT_H

#include <string>
#include <fstream>

#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "../lib/BasicTypes.h" //For ulong, uint
#include "../lib/StrLib.h"
#include "PDBSetup.h" //For atoms class.
#include "EnergyTypes.h"

struct Histogram : OutputableBase
{

  Histogram(OutputVars & v);

  ~Histogram();

  virtual void Sample(const ulong step);

  //No additional init.
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void DoOutput(const ulong step);

private:
  void PrintKindHist(const uint b, const uint k);

  std::string GetFName(std::string const& histName,
                       std::string const& histNum,
                       std::string const& histLetter,
                       const uint box, const uint totKinds);

  //Indices 1: boxes 2: kinds 3: count bins up to N_total
  uint ** molCount[BOXES_WITH_U_NB];
  uint * total;

  std::ofstream * outF [BOXES_WITH_U_NB];
  std::string * name [BOXES_WITH_U_NB];
};

#endif /*HIST_OUTPUT_H*/
