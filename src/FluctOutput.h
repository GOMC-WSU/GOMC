#ifndef FLUCT_OUTPUT_H
#define FLUCT_OUTPUT_H

#include <string>
#include <fstream>

#include "../lib/BasicTypes.h" //For ulong, uint
#include "EnergyTypes.h" //For energies.
#include "MoveSettings.h" //For move settings/state
#include "MoleculeKind.h" //For kind names
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "OutConst.h"
#include "OutputAbstracts.h"
#include "System.h"
#include "StaticVals.h"
#include "PDBSetup.h" //For atoms class.
#include "BoxDimensions.h" //For "boxes with volume."

#include <limits> //for std::numeric_limits

struct FluctuationTracker
{
   FluctuationTracker(): dblSrc(NULL) {}
   
   ~FluctuationTracker() 
   { 
      if (outF.is_open())
      {
	 outF.close();
      }
      if (dblSrc != NULL)
      {
	 delete[] dblSrc;
      }
      if (uintSrc != NULL)
      {
	 delete[] uintSrc;
      }
   }
   
   //Initializes name, and enable
   void Init(const bool en, std::string const& var,
             std::string const& uniqueName, const uint bTot = BOX_TOTAL);

   //Set one of the pointers to the fluctuating values we're tracking
   void SetRef(double * loc, const uint b) 
   {
      dblSrc[b] = loc;
      uintSrc[b] = NULL; 
      outF << std::setprecision(std::numeric_limits<double>::digits10+2) << std::setw(25);
   }
   void SetRef(uint * loc, const uint b) 
   { uintSrc[b] = loc; dblSrc[b] = NULL; }

   void Write(const ulong step, const bool firstPrint)
   { 
      first = firstPrint;
      if (enable)
	 DoWrite(step);
   }

 private:
   
   std::string GetFName(std::string const& var, std::string const& uniqueName);
   
   void DoWrite(const ulong step);

   
   bool first;
   std::ofstream outF;
   std::string name, varName;
   uint ** uintSrc, tot;
   double ** dblSrc;
   bool enable;
};

struct Fluctuations : OutputableBase
{
   Fluctuations(OutputVars & v){ this->var = &v; }

   //Fluctuations does not need to sample, so does nothing.
   virtual void Sample(const ulong step) {}

   //No additional init.
   virtual void Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output);
   
   virtual void DoOutput(const ulong step);
  
 private:   

   void InitVals(config_setup::EventSettings const& event)
   {
      stepsPerOut = event.frequency;
      enableOut = event.enable;
   }
   
   void AllocFlucts(void)
   {
      numKindFlucts = out::TOTAL_K * var->numKinds;
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
      //we don't have mole fraction with only one kind
      if (var->numKinds == 1)
         numKindFlucts--;
#endif
      totalFlucts = out::TOTAL_SINGLE + numKindFlucts;
      flucts= new FluctuationTracker[totalFlucts];
   }
   
   void InitWatchSingle(config_setup::TrackedVars const& tracked);

   void InitWatchMulti(config_setup::TrackedVars const& tracked);
   
   FluctuationTracker * flucts;
   uint numKindFlucts, totalFlucts;
};

#endif /*FLUCT_OUTPUT_H*/
