/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BlockOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "StrLib.h"

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

#include <iostream> //for endl;

#define OUTPUTWIDTH 16

void BlockAverage::Init(std::ofstream* file0,
                        std::ofstream* file1,
                        const bool en,
                        const real scale,
                        std::string const& var,
                        const uint bTot)
{
  outBlock0 = file0;
  outBlock1 = file1;
  tot = bTot;
  block = new real[tot];
  uintSrc = new uint *[tot];
  dblSrc = new real *[tot];
  enable = en;
  scl = scale;
  if (enable) {
    Zero();
    for (uint b = 0; b < tot; ++b) {
      uintSrc[b] = NULL;
      dblSrc[b] = NULL;
    }
    printTitle(var, bTot);
  }
}

void BlockAverage::Sum(void)
{
  if (enable && uintSrc[0] != NULL)
    for (uint b = 0; b < tot; ++b)
      block[b] += (real)(*uintSrc[b]) * scl;
  else if (enable)
    for (uint b = 0; b < tot; ++b)
      block[b] += *dblSrc[b] * scl;
}

void BlockAverage::DoWrite(const ulong step, uint precision)
{
  if (tot >= 1) {
    if (outBlock0->is_open()) {
      if(abs(block[0]) > 1e99) {
        (*outBlock0) << right << std::scientific  << std::setprecision(precision - 1) <<
                     std::setw(OUTPUTWIDTH);
        (*outBlock0) << block[0];
      } else {
        (*outBlock0) << right << std::scientific  << std::setprecision(precision) <<
                     std::setw(OUTPUTWIDTH);
        (*outBlock0) << block[0];
      }
    } else
      std::cerr << "Unable to write to Box_0 output file" << std::endl;
  }
  if (tot >= 2) {
    if (outBlock1->is_open()) {
      if(abs(block[0]) > 1e99) {
        (*outBlock1) << right << std::scientific  << std::setprecision(precision - 1) <<
                     std::setw(OUTPUTWIDTH);
        (*outBlock1) << block[1];
      } else {
        (*outBlock1) << right << std::scientific  << std::setprecision(precision) <<
                     std::setw(OUTPUTWIDTH);
        (*outBlock1) << block[1];
      }
    } else
      std::cerr << "Unable to write to Box_1 output file" << std::endl;
  }
  Zero();
}

void BlockAverages::Init(pdb_setup::Atoms const& atoms,
                         config_setup::Output const& output)
{
  std::string name = "Blk_" + uniqueName + "_BOX_0.dat";
  outBlock0.open(name.c_str(), std::ofstream::out);
  if(BOXES_WITH_U_NB >= 2) {
    name = "Blk_" + uniqueName + "_BOX_1.dat";
    outBlock1.open(name.c_str(), std::ofstream::out);
  }
  InitVals(output.statistics.settings.block);
  AllocBlocks();
  InitWatchSingle(output.statistics.vars);
  InitWatchMulti(output.statistics.vars);
  outBlock0 << std::endl;
  if(outBlock1.is_open())
    outBlock1 << std::endl;
}

void BlockAverages::AllocBlocks(void)
{
  numKindBlocks = out::TOTAL_K * var->numKinds;
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  //we don't have mole fraction and mol density with only one kind
  if (var->numKinds == 1)
    numKindBlocks = 0;
#endif
  totalBlocks = out::TOTAL_SINGLE + numKindBlocks;
  blocks = new BlockAverage[totalBlocks];
}

void BlockAverages::Sample(const ulong step)
{
  for (uint v = 0; v < totalBlocks; ++v)
    blocks[v].Sum();
}

void BlockAverages::DoOutput(const ulong step)
{
  ulong nextStep = step + 1;
  outBlock0 << left << std::scientific << std::setw(OUTPUTWIDTH) << nextStep;
  outBlock1 << left << std::scientific << std::setw(OUTPUTWIDTH) << nextStep;
  for (uint v = 0; v < totalBlocks; ++v) {
    if(v < out::TOTAL_SINGLE)
      blocks[v].Write(nextStep, firstPrint, 8);
    else
      blocks[v].Write(nextStep, firstPrint, 8);
  }
  outBlock0 << std::endl;
  if(outBlock1.is_open())
    outBlock1 << std::endl;
}

void BlockAverages::InitWatchSingle(config_setup::TrackedVars const& tracked)
{
  outBlock0 << left << std::scientific << std::setw(OUTPUTWIDTH) << "#STEPS";
  if(outBlock1.is_open())
    outBlock1 << left << std::scientific << std::setw(OUTPUTWIDTH) << "#STEPS";
  //Note: The order of Init should be same as order of SetRef
  blocks[out::ENERGY_TOTAL_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_TOTAL, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTER_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_INTER, BOXES_WITH_U_NB);
  blocks[out::ENERGY_TC_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_TC, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTRA_B_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_INTRA_B, BOXES_WITH_U_NB);
  blocks[out::ENERGY_INTRA_NB_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_INTRA_NB, BOXES_WITH_U_NB);
  blocks[out::ENERGY_ELECT_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_ELECT, BOXES_WITH_U_NB);
  blocks[out::ENERGY_REAL_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_REAL, BOXES_WITH_U_NB);
  blocks[out::ENERGY_RECIP_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::ENERGY_RECIP, BOXES_WITH_U_NB);
  blocks[out::VIRIAL_TOTAL_IDX].Init(&outBlock0, &outBlock1, tracked.pressure.block, invSteps, out::VIRIAL_TOTAL, BOXES_WITH_U_NB);
  blocks[out::PRESSURE_IDX].Init(&outBlock0, &outBlock1, tracked.pressure.block, invSteps, out::PRESSURE, BOXES_WITH_U_NB);
  blocks[out::MOL_NUM_IDX].Init(&outBlock0, &outBlock1, tracked.molNum.block, invSteps, out::MOL_NUM, BOXES_WITH_U_NB);
  blocks[out::DENSITY_IDX].Init(&outBlock0, &outBlock1, tracked.density.block, invSteps, out::DENSITY, BOXES_WITH_U_NB);
  blocks[out::SURF_TENSION_IDX].Init(&outBlock0, &outBlock1, tracked.surfaceTension.block, invSteps, out::SURF_TENSION, BOXES_WITH_U_NB);
#if ENSEMBLE == GEMC
  blocks[out::VOLUME_IDX].Init(&outBlock0, &outBlock1, tracked.volume.block, invSteps, out::VOLUME, BOXES_WITH_U_NB);
  blocks[out::HEAT_OF_VAP_IDX].Init(&outBlock0, &outBlock1, tracked.energy.block, invSteps, out::HEAT_OF_VAP, BOXES_WITH_U_NB);
#endif
#if ENSEMBLE == NPT
  blocks[out::VOLUME_IDX].Init(&outBlock0, &outBlock1, tracked.volume.block, invSteps, out::VOLUME, BOXES_WITH_U_NB);
#endif

  //Note: The order of Init should be same as order of Init
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    blocks[out::ENERGY_TOTAL_IDX].SetRef(&var->energyRef[b].total, b);
    blocks[out::ENERGY_INTRA_B_IDX].SetRef(&var->energyRef[b].intraBond, b);
    blocks[out::ENERGY_INTER_IDX].SetRef(&var->energyRef[b].inter, b);
    blocks[out::ENERGY_TC_IDX].SetRef(&var->energyRef[b].tc, b);
    blocks[out::ENERGY_INTRA_NB_IDX].SetRef(&var->energyRef[b].intraNonbond, b);
    blocks[out::ENERGY_ELECT_IDX].SetRef(&var->energyRef[b].totalElect, b);
    blocks[out::ENERGY_REAL_IDX].SetRef(&var->energyRef[b].real_en, b);
    blocks[out::ENERGY_RECIP_IDX].SetRef(&var->energyRef[b].recip, b);
    blocks[out::VIRIAL_TOTAL_IDX].SetRef(&var->virialRef[b].total, b);
    blocks[out::PRESSURE_IDX].SetRef(&var->pressure[b], b);
    blocks[out::MOL_NUM_IDX].SetRef(&var->numByBox[b], b);
    blocks[out::DENSITY_IDX].SetRef(&var->densityTot[b], b);
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

void BlockAverages::InitWatchMulti(config_setup::TrackedVars const& tracked)
{
  using namespace pdb_entry::atom::field;
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  uint start = out::TOTAL_SINGLE;
  //Var is molecule kind name plus the prepend related output info kind.
  std::string name;
  for (uint k = 0; k < var->numKinds; ++k) {
    uint bkStart = start + k;
    //Copy each char of the name string.
    std::string trimKindName = var->kindsRef[k].name;
    //If more than one kind, output mol fractions.
    if (var->numKinds > 1) {
      name = out::MOL_FRACTION + "_" + trimKindName;
      blocks[bkStart + out::MOL_FRACTION_IDX * var->numKinds].Init
      (&outBlock0, &outBlock1, tracked.molNum.block, invSteps, name, BOXES_WITH_U_NB);
    }
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      uint kArrIdx = b * var->numKinds + k;
      if (var->numKinds > 1) {
        blocks[bkStart + out::MOL_FRACTION_IDX * var->numKinds].SetRef
        (&var->molFractionByKindBox[kArrIdx], b);
      }
    }
  }
  //I cannot put this in previous loop because the Title will be printed
  // as soon as it is initialized.
  for (uint k = 0; k < var->numKinds; ++k) {
    uint bkStart = start + k;
    //Copy each char of the name string.
    std::string trimKindName = var->kindsRef[k].name;
    //If more than one kind, output mol fractions.
    if (var->numKinds > 1) {
      //Init mol density
      name = out::MOL_DENSITY + "_" + trimKindName;
      blocks[bkStart + out::MOL_DENSITY_IDX * var->numKinds].Init
      (&outBlock0, &outBlock1, tracked.molNum.block, invSteps, name, BOXES_WITH_U_NB);
    }
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      uint kArrIdx = b * var->numKinds + k;
      if (var->numKinds > 1) {
        blocks[bkStart + out::MOL_DENSITY_IDX * var->numKinds].SetRef
        (&var->densityByKindBox[kArrIdx], b);
      }
    }
  }
#endif
}

void BlockAverage::printTitle(std::string output, uint boxes)
{
  if(tot >= 1) {
    if((*outBlock0).is_open()) {
      (*outBlock0) << left << std::scientific << std::setw(OUTPUTWIDTH) << output;
    } else {
      std::cerr << "Unable to write to Block_0 output file!" << std::endl;
    }
  }
  if(tot >= 2) {
    if((*outBlock1).is_open()) {
      (*outBlock1) << left << std::scientific << std::setw(OUTPUTWIDTH) << output;
    } else {
      std::cerr << "Unable to write to Block_1 output file!" << std::endl;
    }
  }
}
