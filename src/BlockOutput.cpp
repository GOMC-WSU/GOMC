#include "BlockOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "../lib/StrLib.h"

#if ENSEMBLE == GEMC
#include "MoveConst.h" //For box constants, if we're calculating Hv
#endif

#include <iostream> //for endl;

void BlockAverage::Init(const bool en, const double scale,
			std::string const& var, std::string const& uniqueName,
			const uint bTot)
{
   tot = bTot;
   block = new double[tot];
   uintSrc = new uint *[tot];
   dblSrc = new double *[tot];
   enable = en;
   scl = scale;
   if (enable)
   {
      name = GetFName(var, uniqueName);
      outF.open(name.c_str(), std::ofstream::out);
      outF.setf(std::ios_base::left, std::ios_base::adjustfield);
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
   if (outF.is_open())
   {
      outF << std::setw(11) << step;
      for (uint b = 0; b < tot; ++b)
      {
	 outF <<  " ";
	 if (dblSrc[b] != NULL)
	    outF << std::setw(25);
	 outF << block[b];
      }
      outF << std::endl;
   }
   else
      std::cerr << "Unable to write to file \"" <<  name << "\" " 
		<< varName << std::endl;
   Zero();
}
   
std::string BlockAverage::GetFName(std::string const& base, 
				   std::string const& uniqueName)
{ 
   std::string fName = "Blk_";
   varName = base;
   fName += base;
   fName += "_";
   fName += uniqueName;
   fName += ".dat";
   return fName;
}

void BlockAverages::Init(pdb_setup::Atoms const& atoms,
                         config_setup::Output const& output)
{
   InitVals(output.statistics.settings.block);
   AllocBlocks();
   InitWatchSingle(output.statistics.vars);
   InitWatchMulti(output.statistics.vars);
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
   for (uint v = 0; v < totalBlocks; ++v)
      blocks[v].Write(nextStep, firstPrint);
}

void BlockAverages::InitWatchSingle(config_setup::TrackedVars const& tracked)
{
   blocks[out::ENERGY_TOTAL_IDX].Init(tracked.energy.block, invSteps, 
				      out::ENERGY_TOTAL, uniqueName, 
				      BOXES_WITH_U_NB);

   //Only output energy categories if specifically requested...
#ifdef EN_SUBCAT_OUT
   blocks[out::ENERGY_INTER_IDX].Init(tracked.energy.block, invSteps, 
				      out::ENERGY_INTER, uniqueName,
				      BOXES_WITH_U_NB);
   blocks[out::ENERGY_TC_IDX].Init(tracked.energy.block, invSteps, 
				   out::ENERGY_TC, uniqueName,
                                   BOXES_WITH_U_NB);
   blocks[out::ENERGY_INTRA_B_IDX].Init(tracked.energy.block, invSteps, 
					out::ENERGY_INTRA_B, uniqueName,
					BOXES_WITH_U_NB);
   blocks[out::ENERGY_INTRA_NB_IDX].Init(tracked.energy.block, invSteps, 
					 out::ENERGY_INTRA_NB, uniqueName,
                                         BOXES_WITH_U_NB);
#endif
   blocks[out::VIRIAL_TOTAL_IDX].Init(tracked.pressure.block, invSteps, 
				      out::VIRIAL_TOTAL, uniqueName,
                                      BOXES_WITH_U_NB);
#ifdef VIR_SUBCAT_OUT
   blocks[out::VIRIAL_INTER_IDX].Init(tracked.pressure.block, invSteps, 
				      out::VIRIAL_INTER, uniqueName,
                                      BOXES_WITH_U_NB);
   blocks[out::VIRIAL_TC_IDX].Init(tracked.pressure.block, invSteps, 
				   out::VIRIAL_TC, uniqueName,
                                   BOXES_WITH_U_NB);
#endif

   blocks[out::PRESSURE_IDX].Init(tracked.pressure.block, invSteps, 
				  out::PRESSURE, uniqueName,
                                  BOXES_WITH_U_NB);
#if ENSEMBLE == GEMC
   blocks[out::VOLUME_IDX].Init(tracked.volume.block, invSteps, 
				out::VOLUME, uniqueName);
   blocks[out::HEAT_OF_VAP_IDX].Init(tracked.energy.block, invSteps,
				     out::HEAT_OF_VAP, uniqueName, 1);

   blocks[out::HEAT_OF_VAP_IDX].SetRef(&var->heatOfVap, 0);
#endif

   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      blocks[out::ENERGY_TOTAL_IDX].SetRef(&var->energyRef[b].total, b);
#ifdef EN_SUBCAT_OUT
      blocks[out::ENERGY_INTRA_B_IDX].SetRef(&var->energyRef[b].intraBond, b);
      blocks[out::ENERGY_INTER_IDX].SetRef(&var->energyRef[b].inter, b);
      blocks[out::ENERGY_TC_IDX].SetRef(&var->energyRef[b].tc, b);
      blocks[out::ENERGY_INTRA_NB_IDX].SetRef
	 (&var->energyRef[b].intraNonbond, b);
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
	 (tracked.molNum.block, invSteps, name, uniqueName, BOXES_WITH_U_NB);
      name = out::DENSITY + "_" + trimKindName;
      blocks[bkStart + out::DENSITY_IDX*var->numKinds].Init
	 (tracked.density.block, invSteps, name, uniqueName, BOXES_WITH_U_NB);
      //If more than one kind, output mol fractions.
      if (var->numKinds > 1)
      {
	 name = out::MOL_FRACTION + "_" + trimKindName;
	 blocks[bkStart + out::MOL_FRACTION_IDX*var->numKinds].Init
	    (tracked.molNum.block, invSteps, name, uniqueName, BOXES_WITH_U_NB);
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
