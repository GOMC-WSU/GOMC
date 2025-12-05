/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "BlockOutput.h"

#include "OutConst.h"
#include "PDBConst.h"
#include "StrLib.h"

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

#include <iostream> //for endl;

#define OUTPUTWIDTH 16

void BlockAverage::Init(std::ofstream *file0, std::ofstream *file1,
                        const bool en, const double scale,
                        const double firstPrint_scale, bool &firstPrint,
                        std::string const &var, const uint bTot) {
  outBlock0 = file0;
  outBlock1 = file1;
  tot = bTot;
  block = new double[tot];
  uintSrc = new uint *[tot];
  dblSrc = new double *[tot];
  enable = en;
  scl = scale;
  fp_scl = firstPrint_scale;
  first = &firstPrint;
  if (enable) {
    Zero();
    for (uint b = 0; b < tot; ++b) {
      uintSrc[b] = NULL;
      dblSrc[b] = NULL;
    }
    printTitle(var);
  }
}

void BlockAverage::Sum(void) {
  if (enable && uintSrc[0] != NULL) {
    for (uint b = 0; b < tot; ++b) {
      // We use a pointer because you can't pass
      // references in the constructor since
      // initializing an array of objects
      // requires the default constructor.
      // This could be fixed by using vectors..
      if (*first)
        block[b] += (double)(*uintSrc[b]) * fp_scl;
      else
        block[b] += (double)(*uintSrc[b]) * scl;
    }
  } else if (enable) {
    for (uint b = 0; b < tot; ++b) {
      // We use a pointer because you can't pass
      // references in the constructor since
      // initializing an array of objects
      // requires the default constructor.
      // This could be fixed by using vectors..
      if (*first)
        block[b] += *dblSrc[b] * fp_scl;
      else
        block[b] += *dblSrc[b] * scl;
    }
  }
}

void BlockAverage::DoWrite(uint precision) {
  if (tot >= 1) {
    if (outBlock0->is_open()) {
      if (std::abs(block[0]) > 1e99) {
        (*outBlock0) << std::right << std::scientific
                     << std::setprecision(precision - 1)
                     << std::setw(OUTPUTWIDTH);
        (*outBlock0) << block[0];
      } else {
        (*outBlock0) << std::right << std::scientific
                     << std::setprecision(precision) << std::setw(OUTPUTWIDTH);
        (*outBlock0) << block[0];
      }
    } else
      std::cerr << "Unable to write to Box_0 output file" << std::endl;
  }
  if (tot >= 2) {
    if (outBlock1->is_open()) {
      if (std::abs(block[0]) > 1e99) {
        (*outBlock1) << std::right << std::scientific
                     << std::setprecision(precision - 1)
                     << std::setw(OUTPUTWIDTH);
        (*outBlock1) << block[1];
      } else {
        (*outBlock1) << std::right << std::scientific
                     << std::setprecision(precision) << std::setw(OUTPUTWIDTH);
        (*outBlock1) << block[1];
      }
    } else
      std::cerr << "Unable to write to Box_1 output file" << std::endl;
  }
  Zero();
}

void BlockAverages::Init(pdb_setup::Atoms const &atoms,
                         config_setup::Output const &output) {
#if GOMC_LIB_MPI
  std::string name =
      pathToReplicaOutputDirectory + "Blk_" + uniqueName + "_BOX_0.dat";
#else
  std::string name = "Blk_" + uniqueName + "_BOX_0.dat";
#endif
  outBlock0.open(name.c_str(), std::ofstream::out);
  if (BOXES_WITH_U_NB >= 2) {
#if GOMC_LIB_MPI
    name = pathToReplicaOutputDirectory + "Blk_" + uniqueName + "_BOX_1.dat";
#else
    name = "Blk_" + uniqueName + "_BOX_1.dat";
#endif
    outBlock1.open(name.c_str(), std::ofstream::out);
  }
  InitVals(output.statistics.settings.block);
  AllocBlocks();
  InitWatchSingle(output.statistics.vars);
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  InitWatchMulti(output.statistics.vars);
#endif
  outBlock0 << std::endl;
  if (outBlock1.is_open())
    outBlock1 << std::endl;
}

void BlockAverages::AllocBlocks(void) {
  numKindBlocks = out::TOTAL_K * var->numKinds;
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  // we don't have mole fraction and mol density with only one kind
  if (var->numKinds == 1)
    numKindBlocks = 0;
#endif
  totalBlocks = out::TOTAL_SINGLE + numKindBlocks;
  blocks = new BlockAverage[totalBlocks];
}

void BlockAverages::Sample(const ulong step) {
  for (uint v = 0; v < totalBlocks; ++v)
    blocks[v].Sum();
}

void BlockAverages::DoOutput(const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::BLK_OUTPUT);
  ulong nextStep = step + 1;
  outBlock0 << std::left << std::scientific << std::setw(OUTPUTWIDTH)
            << nextStep;
  outBlock1 << std::left << std::scientific << std::setw(OUTPUTWIDTH)
            << nextStep;
  for (uint v = 0; v < totalBlocks; ++v) {
    if (v < out::TOTAL_SINGLE)
      blocks[v].Write(8);
    else
      blocks[v].Write(8);
  }
  outBlock0 << std::endl;
  if (outBlock1.is_open())
    outBlock1 << std::endl;
  GOMC_EVENT_STOP(1, GomcProfileEvent::BLK_OUTPUT);
}

void BlockAverages::DoOutputRestart(const ulong step) {}

void BlockAverages::InitWatchSingle(config_setup::TrackedVars const &tracked) {
  outBlock0 << std::left << std::scientific << std::setw(OUTPUTWIDTH)
            << "#STEPS";
  if (outBlock1.is_open())
    outBlock1 << std::left << std::scientific << std::setw(OUTPUTWIDTH)
              << "#STEPS";
  // Note: The order of Init should be same as order of SetRef
  blocks[out::ENERGY_TOTAL_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_TOTAL, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTER_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_INTER, BOXES_WITH_U_NB);
  blocks[out::ENERGY_LRC_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block,
                                   invSteps, firstInvSteps, firstPrint,
                                   out::ENERGY_LRC, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTRA_B_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_INTRA_B, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTRA_NB_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_INTRA_NB, BOXES_WITH_U_NB);
  blocks[out::ENERGY_ELECT_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_ELECT, BOXES_WITH_U_NB);
  blocks[out::ENERGY_REAL_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_REAL, BOXES_WITH_U_NB);
  blocks[out::ENERGY_RECIP_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::ENERGY_RECIP, BOXES_WITH_U_NB);
  blocks[out::VIRIAL_TOTAL_IDX].Init(
      &outBlock0, &outBlock1, tracked.pressure.block, invSteps, firstInvSteps,
      firstPrint, out::VIRIAL_TOTAL, BOXES_WITH_U_NB);
  blocks[out::PRESSURE_IDX].Init(&outBlock0, &outBlock1, tracked.pressure.block,
                                 invSteps, firstInvSteps, firstPrint,
                                 out::PRESSURE, BOXES_WITH_U_NB);
  blocks[out::MOL_NUM_IDX].Init(&outBlock0, &outBlock1, tracked.molNum.block,
                                invSteps, firstInvSteps, firstPrint,
                                out::MOL_NUM, BOXES_WITH_U_NB);
  blocks[out::DENSITY_IDX].Init(&outBlock0, &outBlock1, tracked.density.block,
                                invSteps, firstInvSteps, firstPrint,
                                out::DENSITY, BOXES_WITH_U_NB);
  blocks[out::COMPRESSIBILITY_IDX].Init(
      &outBlock0, &outBlock1, tracked.pressure.block, invSteps, firstInvSteps,
      firstPrint, out::COMPRESSIBILITY, BOXES_WITH_U_NB);
  blocks[out::ENTHALPY_IDX].Init(&outBlock0, &outBlock1, tracked.pressure.block,
                                 invSteps, firstInvSteps, firstPrint,
                                 out::ENTHALPY, BOXES_WITH_U_NB);
  blocks[out::SURF_TENSION_IDX].Init(
      &outBlock0, &outBlock1, tracked.surfaceTension.block, invSteps,
      firstInvSteps, firstPrint, out::SURF_TENSION, BOXES_WITH_U_NB);
#if ENSEMBLE == GEMC
  blocks[out::VOLUME_IDX].Init(&outBlock0, &outBlock1, tracked.volume.block,
                               invSteps, firstInvSteps, firstPrint, out::VOLUME,
                               BOXES_WITH_U_NB);
  blocks[out::HEAT_OF_VAP_IDX].Init(
      &outBlock0, &outBlock1, tracked.energy.block, invSteps, firstInvSteps,
      firstPrint, out::HEAT_OF_VAP, BOXES_WITH_U_NB);
#endif
#if ENSEMBLE == NPT
  blocks[out::VOLUME_IDX].Init(&outBlock0, &outBlock1, tracked.volume.block,
                               invSteps, firstInvSteps, firstPrint, out::VOLUME,
                               BOXES_WITH_U_NB);
#endif

  // Note: The order of Init should be same as order of Init
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    blocks[out::ENERGY_TOTAL_IDX].SetRef(&var->energyRef[b].total, b);
    blocks[out::ENERGY_INTRA_B_IDX].SetRef(&var->energyRef[b].intraBond, b);
    blocks[out::ENERGY_INTER_IDX].SetRef(&var->energyRef[b].inter, b);
    blocks[out::ENERGY_LRC_IDX].SetRef(&var->energyRef[b].tailCorrection, b);
    blocks[out::ENERGY_INTRA_NB_IDX].SetRef(&var->energyRef[b].intraNonbond, b);
    blocks[out::ENERGY_ELECT_IDX].SetRef(&var->energyRef[b].totalElect, b);
    blocks[out::ENERGY_REAL_IDX].SetRef(&var->energyRef[b].real, b);
    blocks[out::ENERGY_RECIP_IDX].SetRef(&var->energyRef[b].recip, b);
    blocks[out::VIRIAL_TOTAL_IDX].SetRef(&var->virialRef[b].total, b);
    blocks[out::PRESSURE_IDX].SetRef(&var->pressure[b], b);
    blocks[out::MOL_NUM_IDX].SetRef(&var->numByBox[b], b);
    blocks[out::DENSITY_IDX].SetRef(&var->densityTot[b], b);
    blocks[out::COMPRESSIBILITY_IDX].SetRef(&var->compressibility[b], b);
    blocks[out::ENTHALPY_IDX].SetRef(&var->enthalpy[b], b);
    blocks[out::SURF_TENSION_IDX].SetRef(&var->surfaceTens[b], b);
#if ENSEMBLE == GEMC
    blocks[out::VOLUME_IDX].SetRef(&var->volumeRef[b], b);
    blocks[out::HEAT_OF_VAP_IDX].SetRef(&var->heatOfVap, b);
#endif
#if ENSEMBLE == NPT
    blocks[out::VOLUME_IDX].SetRef(&var->volumeRef[b], b);
#endif
  }
}

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
void BlockAverages::InitWatchMulti(config_setup::TrackedVars const &tracked) {
  using namespace pdb_entry::atom::field;
  uint start = out::TOTAL_SINGLE;
  // Var is molecule kind name plus the prepend related output info kind.
  std::string name;
  for (uint k = 0; k < var->numKinds; ++k) {
    uint bkStart = start + k;
    // Copy each char of the name string.
    std::string trimKindName = var->kindsRef[k].name;
    // If more than one kind, output mol fractions.
    if (var->numKinds > 1) {
      name = out::MOL_FRACTION + "_" + trimKindName;
      blocks[bkStart + out::MOL_FRACTION_IDX * var->numKinds].Init(
          &outBlock0, &outBlock1, tracked.molNum.block, invSteps, firstInvSteps,
          firstPrint, name, BOXES_WITH_U_NB);
    }
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      uint kArrIdx = b * var->numKinds + k;
      if (var->numKinds > 1) {
        blocks[bkStart + out::MOL_FRACTION_IDX * var->numKinds].SetRef(
            &var->molFractionByKindBox[kArrIdx], b);
      }
    }
  }
  // I cannot put this in previous loop because the Title will be printed
  // as soon as it is initialized.
  for (uint k = 0; k < var->numKinds; ++k) {
    uint bkStart = start + k;
    // Copy each char of the name string.
    std::string trimKindName = var->kindsRef[k].name;
    // If more than one kind, output mol fractions.
    if (var->numKinds > 1) {
      // Init mol density
      name = out::MOL_DENSITY + "_" + trimKindName;
      blocks[bkStart + out::MOL_DENSITY_IDX * var->numKinds].Init(
          &outBlock0, &outBlock1, tracked.molNum.block, invSteps, firstInvSteps,
          firstPrint, name, BOXES_WITH_U_NB);
    }
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      uint kArrIdx = b * var->numKinds + k;
      if (var->numKinds > 1) {
        blocks[bkStart + out::MOL_DENSITY_IDX * var->numKinds].SetRef(
            &var->densityByKindBox[kArrIdx], b);
      }
    }
  }
}
#endif

void BlockAverage::printTitle(std::string output) {
  if (tot >= 1) {
    if ((*outBlock0).is_open()) {
      (*outBlock0) << std::left << std::scientific << std::setw(OUTPUTWIDTH)
                   << output;
    } else {
      std::cerr << "Unable to write to Block_0 output file!" << std::endl;
    }
  }
  if (tot >= 2) {
    if ((*outBlock1).is_open()) {
      (*outBlock1) << std::left << std::scientific << std::setw(OUTPUTWIDTH)
                   << output;
    } else {
      std::cerr << "Unable to write to Block_1 output file!" << std::endl;
    }
  }
}