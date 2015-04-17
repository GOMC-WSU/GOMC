/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FluctOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "../lib/StrLib.h"

void FluctuationTracker::Init(const bool en, std::string const& var,
                              std::string const& uniqueName, const uint bTot)
{
   tot = bTot;
   uintSrc = new uint *[tot];
   dblSrc = new double *[tot];
   first = enable = en;
   if (enable)
      name = GetFName(var, uniqueName);
}

void FluctuationTracker::DoWrite(const ulong step)
{
   outF.open(name.c_str(), (first?std::ofstream::out:std::ofstream::app));
   if (outF.is_open())
   {
      outF << step;
      for (uint b = 0; b < tot; ++b)
         if (uintSrc[0] != NULL)
            outF << " " << *uintSrc[b];
         else
            outF << " " << *dblSrc[b];
      outF << std::endl;
      outF.close();
   }
   else
      std::cerr << "Unable to write to file \"" <<  name << "\" " 
		<< varName << std::endl;
}
   
std::string FluctuationTracker::GetFName(std::string const& base, 
                                         std::string const& uniqueName)
{ 
   std::string fName = "Fluct_";
   varName = base;
   fName += base;
   fName += "_";
   fName += uniqueName;
   fName += ".dat";
   return fName;
}


void Fluctuations::Init(pdb_setup::Atoms const& atoms,
                        config_setup::Output const& output)
{
   InitVals(output.statistics.settings.fluct);
   AllocFlucts();
   InitWatchSingle(output.statistics.vars);
   InitWatchMulti(output.statistics.vars);
}

void Fluctuations::DoOutput(const ulong step)
{
   ulong nextStep = step+1;
   for (uint v = 0; v < totalFlucts; ++v)
      flucts[v].Write(nextStep, firstPrint);
}

void Fluctuations::InitWatchSingle(config_setup::TrackedVars const& tracked)
{
   flucts[out::ENERGY_TOTAL_IDX].Init(tracked.energy.fluct, out::ENERGY_TOTAL,
                                      uniqueName);
   flucts[out::ENERGY_INTER_IDX].Init(tracked.energy.fluct,
				      out::ENERGY_INTER, uniqueName,
				      BOXES_WITH_U_NB);
   flucts[out::ENERGY_TC_IDX].Init(tracked.energy.fluct,
				   out::ENERGY_TC, uniqueName,
				   BOXES_WITH_U_NB);
   flucts[out::ENERGY_INTRA_B_IDX].Init(tracked.energy.fluct, 
					out::ENERGY_INTRA_B, uniqueName);
   flucts[out::ENERGY_INTRA_NB_IDX].Init(tracked.energy.fluct, 
					 out::ENERGY_INTRA_NB, uniqueName,
					 BOXES_WITH_U_NB);
   flucts[out::VIRIAL_TOTAL_IDX].Init(tracked.pressure.fluct, 
				      out::VIRIAL_TOTAL, uniqueName,
					 BOXES_WITH_U_NB);
   flucts[out::VIRIAL_INTER_IDX].Init(tracked.pressure.fluct, 
				      out::VIRIAL_INTER, uniqueName,
					 BOXES_WITH_U_NB);
   flucts[out::VIRIAL_TC_IDX].Init(tracked.pressure.fluct, 
				   out::VIRIAL_TC, uniqueName,
				   BOXES_WITH_U_NB);
   flucts[out::PRESSURE_IDX].Init(tracked.pressure.fluct, 
				  out::PRESSURE, uniqueName,
				  BOXES_WITH_U_NB);
#if ENSEMBLE == GEMC
   flucts[out::VOLUME_IDX].Init(tracked.volume.fluct,
				out::VOLUME, uniqueName);
   flucts[out::HEAT_OF_VAP_IDX].Init(tracked.energy.fluct,
				     out::HEAT_OF_VAP, uniqueName, 1);

   flucts[out::HEAT_OF_VAP_IDX].SetRef(&var->heatOfVap, 0);
#endif
   
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      flucts[out::ENERGY_TOTAL_IDX].SetRef(&var->energyRef[b].total, b);
      flucts[out::ENERGY_INTRA_B_IDX].SetRef(&var->energyRef[b].intraBond, b);
#if ENSEMBLE == GEMC
      flucts[out::VOLUME_IDX].SetRef(&var->volumeRef[b], b);
#endif
   }

   for (uint b = 0; b < BOXES_WITH_U_NB; ++b)
   {
      flucts[out::ENERGY_INTER_IDX].SetRef(&var->energyRef[b].inter, b);
      flucts[out::ENERGY_TC_IDX].SetRef(&var->energyRef[b].tc, b);
      flucts[out::ENERGY_INTRA_NB_IDX].SetRef
	 (&var->energyRef[b].intraNonbond, b);
      flucts[out::VIRIAL_TOTAL_IDX].SetRef(&var->virialRef[b].total, b);
      flucts[out::VIRIAL_INTER_IDX].SetRef(&var->virialRef[b].inter, b);
      flucts[out::VIRIAL_TC_IDX].SetRef(&var->virialRef[b].tc, b);
      flucts[out::PRESSURE_IDX].SetRef(&var->pressure[b], b);
   } 
}

void Fluctuations::InitWatchMulti(config_setup::TrackedVars const& tracked)
{
   using namespace pdb_entry::atom::field;
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
   uint start = out::TOTAL_SINGLE;
   //Var is molecule kind name plus the prepend related output info kind.
   std::string name;
   for (uint k = 0; k < var->numKinds; k++)
   {
      uint fkStart = start + k;
      //Copy each char of the name string.
      std::string trimKindName = var->kindsRef[k].name;
      name = out::MOL_NUM + "_" + trimKindName;
      flucts[fkStart + out::MOL_NUM_IDX*var->numKinds].Init
	 (tracked.molNum.fluct, name, uniqueName);
      name = out::DENSITY + "_" + trimKindName;
      flucts[fkStart + out::DENSITY_IDX*var->numKinds].Init
	 (tracked.density.fluct, name, uniqueName);
      //If more than one kind, output mol fractions.
      if (var->numKinds > 1)
      {
	 name = out::MOL_FRACTION + "_" + trimKindName;
	 flucts[fkStart + out::MOL_FRACTION_IDX*var->numKinds].Init
	    (tracked.molNum.fluct, name, uniqueName);
      }
      for (uint b = 0; b < BOX_TOTAL; ++b)
      {
	 uint kArrIdx = b*var->numKinds+k;
	 flucts[fkStart + out::MOL_NUM_IDX*var->numKinds].SetRef
	    (&var->numByKindBox[kArrIdx], b);
	 flucts[fkStart + out::DENSITY_IDX*var->numKinds].SetRef
            (&var->densityByKindBox[kArrIdx], b);
	 if (var->numKinds > 1)
	    flucts[fkStart + out::MOL_FRACTION_IDX*var->numKinds].SetRef
	       (&var->molFractionByKindBox[kArrIdx], b);
      }
   }
#endif
}

