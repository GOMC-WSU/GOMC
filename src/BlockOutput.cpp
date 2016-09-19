#include "BlockOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "../lib/StrLib.h"

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

#include <iostream> //for endl;

void BlockAverage::Init(std::ofstream* file,
			const bool en, 
			const double scale,
                        std::string const& var,
                        const uint bTot)
{
  outF = file;
  tot = bTot;
  block = new double[tot];
  uintSrc = new uint *[tot];
  dblSrc = new double *[tot];
  enable = en;
  scl = scale;
  if (enable)
  {
    Zero();
    for (uint b = 0; b < tot; ++b)
    {
      uintSrc[b] = NULL;
      dblSrc[b] = NULL;
    }
  }
}

void BlockAverage::Sum(void)
{
  if (enable && uintSrc[0] != NULL)
    for (uint b = 0; b < tot; ++b)
      block[b] += (double)(*uintSrc[b]) * scl;
  else if (enable)
    for (uint b = 0; b < tot; ++b)
      block[b] += *dblSrc[b] * scl;
}

void BlockAverage::DoWrite(const ulong step)
{
  if (outF->is_open())
  {
    for (uint b = 0; b < tot; ++b)
    {
      (*outF) << std::setw(25);
      (*outF) << block[b];
    }
  }
  else
    std::cerr << "Unable to write to file \"" <<  name << "\" "
              << varName << std::endl;
  Zero();
}

void BlockAverages::Init(pdb_setup::Atoms const& atoms,
                         config_setup::Output const& output)
{
  std::string name = "Blk_Average_Energy.dat";
  outF.open(name.c_str(), std::ofstream::out);    
  InitVals(output.statistics.settings.block);
  AllocBlocks();
  InitWatchSingle(output.statistics.vars);
  InitWatchMulti(output.statistics.vars);
  outF << std::endl;
}

void BlockAverages::AllocBlocks(void)
{
  numKindBlocks = out::TOTAL_K * var->numKinds;
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  //we don't have mole fraction with only one kind
  if (var->numKinds == 1)
    numKindBlocks--;
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
  ulong nextStep = step+1;
  outF << std::setw(25) << step;
  for (uint v = 0; v < totalBlocks; ++v)
    blocks[v].Write(nextStep, firstPrint);
  outF << std::endl;
}

void BlockAverages::InitWatchSingle(config_setup::TrackedVars const& tracked)
{
  
  outF << std::setw(25) << "TITLES";
  blocks[out::ENERGY_TOTAL_IDX].Init(&outF, tracked.energy.block, invSteps,
                                     out::ENERGY_TOTAL, BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_TOTAL;

  //Only output energy categories if specifically requested...
#ifdef EN_SUBCAT_OUT
  blocks[out::ENERGY_INTER_IDX].Init(&outF, tracked.energy.block, invSteps,
                                     out::ENERGY_INTER,
                                     BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_INTER;
  blocks[out::ENERGY_TC_IDX].Init(&outF, tracked.energy.block, invSteps,
                                  out::ENERGY_TC,
                                  BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_TC;
  blocks[out::ENERGY_INTRA_B_IDX].Init(&outF, tracked.energy.block, invSteps,
                                       out::ENERGY_INTRA_B,
                                       BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_INTRA_B;
  blocks[out::ENERGY_INTRA_NB_IDX].Init(&outF, tracked.energy.block, invSteps,
                                        out::ENERGY_INTRA_NB,
                                        BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_INTRA_NB;
  blocks[out::ENERGY_ELECT_IDX].Init(&outF, tracked.energy.block, invSteps,
                                     out::ENERGY_ELECT,
                                     BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_ELECT;
  blocks[out::ENERGY_REAL_IDX].Init(&outF, tracked.energy.block, invSteps,
                                    out::ENERGY_REAL,
                                    BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_REAL;
  blocks[out::ENERGY_RECIP_IDX].Init(&outF, tracked.energy.block, invSteps,
                                     out::ENERGY_RECIP,
                                     BOXES_WITH_U_NB);
  outF << std::setw(25) << out::ENERGY_RECIP;
#endif
  blocks[out::VIRIAL_TOTAL_IDX].Init(&outF, tracked.pressure.block, invSteps,
                                     out::VIRIAL_TOTAL,
                                     BOXES_WITH_U_NB);
  outF << std::setw(25) << out::VIRIAL_TOTAL;
#ifdef VIR_SUBCAT_OUT
  blocks[out::VIRIAL_INTER_IDX].Init(&outF, tracked.pressure.block, invSteps,
                                     out::VIRIAL_INTER,
                                     BOXES_WITH_U_NB);
  outF << std::setw(25) << out::VIRIAL_INTER;
  blocks[out::VIRIAL_TC_IDX].Init(&outF, tracked.pressure.block, invSteps,
                                  out::VIRIAL_TC,
                                  BOXES_WITH_U_NB);
  outF << std::setw(25) << out::VIRIAL_TC;
#endif

  blocks[out::PRESSURE_IDX].Init(&outF, tracked.pressure.block, invSteps,
                                 out::PRESSURE,
                                 BOXES_WITH_U_NB);
  outF << std::setw(25) << out::PRESSURE;
#if ENSEMBLE == GEMC
  blocks[out::VOLUME_IDX].Init(&outF, tracked.volume.block, invSteps,
                               out::VOLUME);
  outF << std::setw(25) << out::VOLUME;
  blocks[out::HEAT_OF_VAP_IDX].Init(&outF, tracked.energy.block, invSteps,
                                    out::HEAT_OF_VAP, 1);
  outF << std::setw(25) << out::HEAT_OF_VAP;

  blocks[out::HEAT_OF_VAP_IDX].SetRef(&var->heatOfVap, 0);
#endif

  for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
  {
    blocks[out::ENERGY_TOTAL_IDX].SetRef(&var->energyRef[b].total, b);
#ifdef EN_SUBCAT_OUT
    blocks[out::ENERGY_INTRA_B_IDX].SetRef(&var->energyRef[b].intraBond, b);
    blocks[out::ENERGY_INTER_IDX].SetRef(&var->energyRef[b].inter, b);
    blocks[out::ENERGY_TC_IDX].SetRef(&var->energyRef[b].tc, b);
    blocks[out::ENERGY_INTRA_NB_IDX].SetRef(&var->energyRef[b].intraNonbond, b);
    blocks[out::ENERGY_ELECT_IDX].SetRef(&var->energyRef[b].totalElect, b);
    blocks[out::ENERGY_REAL_IDX].SetRef(&var->energyRef[b].real, b);
    blocks[out::ENERGY_RECIP_IDX].SetRef(&var->energyRef[b].recip, b);
#endif
    blocks[out::VIRIAL_TOTAL_IDX].SetRef(&var->virialRef[b].total, b);
#ifdef VIR_SUBCAT_OUT
    blocks[out::VIRIAL_INTER_IDX].SetRef(&var->virialRef[b].inter, b);
    blocks[out::VIRIAL_TC_IDX].SetRef(&var->virialRef[b].tc, b);
#endif
    blocks[out::PRESSURE_IDX].SetRef(&var->pressure[b], b);
#if ENSEMBLE == GEMC
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
  for (uint k = 0; k < var->numKinds; ++k)
  {
    uint bkStart = start + k;
    //Copy each char of the name string.
    std::string trimKindName = var->kindsRef[k].name;
    name = out::MOL_NUM + "_" + trimKindName;
    blocks[bkStart + out::MOL_NUM_IDX*var->numKinds].Init
      (&outF, tracked.molNum.block, invSteps, name, BOXES_WITH_U_NB);
    outF << std::setw(25) << name;
    name = out::DENSITY + "_" + trimKindName;
    blocks[bkStart + out::DENSITY_IDX*var->numKinds].Init
      (&outF, tracked.density.block, invSteps, name, BOXES_WITH_U_NB);
    outF << std::setw(25) << name;
    //If more than one kind, output mol fractions.
    if (var->numKinds > 1)
    {
      name = out::MOL_FRACTION + "_" + trimKindName;
      blocks[bkStart + out::MOL_FRACTION_IDX*var->numKinds].Init
	(&outF, tracked.molNum.block, invSteps, name, BOXES_WITH_U_NB);
      outF << std::setw(25) << name;
    }    
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
    {
      uint kArrIdx = b*var->numKinds+k;
      blocks[bkStart + out::MOL_NUM_IDX*var->numKinds].SetRef
      (&var->numByKindBox[kArrIdx], b);
      blocks[bkStart + out::DENSITY_IDX*var->numKinds].SetRef
      (&var->densityByKindBox[kArrIdx], b);
      if (var->numKinds > 1)
        blocks[bkStart + out::MOL_FRACTION_IDX*var->numKinds].SetRef
        (&var->molFractionByKindBox[kArrIdx], b);
    }
  }
#endif
}
