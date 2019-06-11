/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BLOCK_OUTPUT_H
#define BLOCK_OUTPUT_H

#include <string>
#include <fstream>

#include "BasicTypes.h" //For ulong, uint
#include "EnergyTypes.h" //For energies.
#include "MoleculeKind.h" //For kind names
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "OutConst.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "StaticVals.h"
#include "PDBSetup.h" //For atoms class.
#include "BoxDimensions.h" //For BOXES_WITH_VOLUME
#include "BoxDimensionsNonOrth.h"

#include <limits> //for std::numeric_limits

class System;

struct BlockAverage {
  BlockAverage(): enable(false), block(NULL), uintSrc(NULL), dblSrc(NULL) {}
  ~BlockAverage()
  {
    if (dblSrc != NULL) {
      delete[] dblSrc;
    }
    if (uintSrc != NULL) {
      delete[] uintSrc;
    }
    if (block != NULL) {
      delete[] block;
    }
  }
  //Initializes name, and enable
  void Init(std::ofstream *file0,
            std::ofstream *file1,
            const bool en,
            const real scl,
            std::string const& var,
            const uint bTot = BOX_TOTAL);

  //Set one of the pointers to the block values we're tracking
  void SetRef(real * loc, const uint b)
  {
    dblSrc[b] = loc;
    uintSrc[b] = NULL;
  }
  void SetRef(uint * loc, const uint b)
  {
    uintSrc[b] = loc;
    dblSrc[b] = NULL;
  }
  void Sum(void);
  void Write(const ulong step, const bool firstPrint, uint precision)
  {
    first = firstPrint;
    if (enable)
      DoWrite(step, precision);
  }

private:
  void Zero(void)
  {
    for (uint b = 0; b < tot; b++)
      block[b] = 0.0;
    samples = 0;
  }
  void DoWrite(const ulong step, uint precision);
  void printTitle(std::string output, uint boxes);

  std::ofstream* outBlock0;
  std::ofstream* outBlock1;
  bool first;
  std::string name, varName;
  uint ** uintSrc, tot;
  real ** dblSrc;
  real * block, scl;
  uint samples;
  bool enable;
};
/**********************************************************************/
struct BlockAverages : OutputableBase {
  BlockAverages(): blocks(NULL) {}

  BlockAverages(OutputVars & v)
  {
    this->var = &v;
    blocks = NULL;
  }

  ~BlockAverages(void)
  {
    if (outBlock0.is_open()) {
      outBlock0.close();
    }
    if (outBlock1.is_open()) {
      outBlock1.close();
    }
    if ( blocks != NULL ) {
      delete[] blocks;
    }
  }
  //No additional init.
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void Sample(const ulong step);
  virtual void DoOutput(const ulong step);

private:
  void InitVals(config_setup::EventSettings const& event)
  {
    stepsPerOut = event.frequency;
    invSteps = 1.0 / stepsPerOut;
    enableOut = event.enable;
  }
  void AllocBlocks(void);
  void InitWatchSingle(config_setup::TrackedVars const& tracked);
  void InitWatchMulti(config_setup::TrackedVars const& tracked);

  std::ofstream outBlock0;
  std::ofstream outBlock1;
  //Block vars
  BlockAverage * blocks;
  uint numKindBlocks, totalBlocks;
  //Intermediate vars.
  uint samplesWrites;
  //Constants
  real invSteps;
};

#endif /*BLOCK_OUTPUT_H*/
