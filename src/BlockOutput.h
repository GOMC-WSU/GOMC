/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef BLOCK_OUTPUT_H
#define BLOCK_OUTPUT_H

#include <fstream>
#include <limits> //for std::numeric_limits
#include <string>

#include "BasicTypes.h"    //For ulong, uint
#include "BoxDimensions.h" //For BOXES_WITH_VOLUME
#include "BoxDimensionsNonOrth.h"
#include "EnergyTypes.h"    //For energies.
#include "MoleculeKind.h"   //For kind names
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "OutConst.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "PDBSetup.h" //For atoms class.
#include "StaticVals.h"

class System;

struct BlockAverage {
  BlockAverage() : uintSrc(NULL), dblSrc(NULL), block(NULL), enable(false) {}
  ~BlockAverage() {
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
  // Initializes name, and enable
  void Init(std::ofstream *file0, std::ofstream *file1, const bool en,
            const double scl, const double firstPrint_scl, bool &firstPrint,
            std::string const &var, const uint bTot = BOX_TOTAL);

  // Set one of the pointers to the block values we're tracking
  void SetRef(double *loc, const uint b) {
    dblSrc[b] = loc;
    uintSrc[b] = NULL;
  }
  void SetRef(uint *loc, const uint b) {
    uintSrc[b] = loc;
    dblSrc[b] = NULL;
  }
  void Sum(void);
  void Write(uint precision) {
    if (enable)
      DoWrite(precision);
  }

private:
  void Zero(void) {
    for (uint b = 0; b < tot; b++)
      block[b] = 0.0;
    samples = 0;
  }
  void DoWrite(uint precision);
  void printTitle(std::string output);

  std::ofstream *outBlock0;
  std::ofstream *outBlock1;
  bool *first;
  std::string name, varName;
  uint **uintSrc, tot;
  double **dblSrc;
  double *block, scl, fp_scl;
  uint samples = 0;
  bool enable;
};
/**********************************************************************/
struct BlockAverages : OutputableBase {
  BlockAverages() : blocks(NULL) {}

  BlockAverages(OutputVars &v) {
    this->var = &v;
    blocks = NULL;
  }

  ~BlockAverages(void) {
    if (outBlock0.is_open()) {
      outBlock0.close();
    }
    if (outBlock1.is_open()) {
      outBlock1.close();
    }
    if (blocks != NULL) {
      delete[] blocks;
    }
  }
  // No additional init.
  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void Sample(const ulong step);
  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  void InitVals(config_setup::EventSettings const &event) {
    stepsPerOut = event.frequency;
    invSteps = 1.0 / stepsPerOut;
    firstInvSteps = invSteps;
    //Handle the case where we are restarting from a checkpoint and the first
    //interval is smaller than expected because we create a checkpoint more
    //often than the Block output frequency.
    if (startStep != 0 && (startStep % stepsPerOut) != 0) {
      ulong diff;
      diff = stepsPerOut - (startStep % stepsPerOut);
      firstInvSteps = 1.0/diff;
    }
    enableOut = event.enable;
  }
  void AllocBlocks(void);
  void InitWatchSingle(config_setup::TrackedVars const &tracked);
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  void InitWatchMulti(config_setup::TrackedVars const &tracked);
#endif
  std::ofstream outBlock0;
  std::ofstream outBlock1;
  // Block vars
  BlockAverage *blocks;
  uint numKindBlocks, totalBlocks;
  // Intermediate vars.
  uint samplesWrites;
  // Constants
  double invSteps;
  double firstInvSteps;
};

#endif /*BLOCK_OUTPUT_H*/
