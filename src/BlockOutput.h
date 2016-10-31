/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BLOCK_OUTPUT_H
#define BLOCK_OUTPUT_H

#include <string>
#include <fstream>

#include "../lib/BasicTypes.h" //For ulong, uint
#include "EnergyTypes.h" //For energies.
#include "MoleculeKind.h" //For kind names
#include "MoleculeLookup.h" //for lookup array (to get density, kind cnts, etc.
#include "OutConst.h"
#include "OutputAbstracts.h"
#include "OutputVars.h"
#include "StaticVals.h"
#include "PDBSetup.h" //For atoms class.
#include "BoxDimensions.h" //For BOXES_WITH_VOLUME

#include <limits> //for std::numeric_limits

class System;

struct BlockAverage
{
   BlockAverage(): enable(false), block(NULL), uintSrc(NULL), dblSrc(NULL) {}

   ~BlockAverage()
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
      if (block != NULL)
      {
	 delete[] block;
      }
   }

   //Initializes name, and enable
   void Init(const bool en, const double scl,
	     std::string const& var, std::string const& uniqueName,
	     const uint bTot = BOX_TOTAL);

   //Set one of the pointers to the block values we're tracking
   void SetRef(double * loc, const uint b)
   {
      dblSrc[b] = loc;
      uintSrc[b] = NULL;
      outF << std::setprecision(std::numeric_limits<double>::digits10+2) << std::setw(25);
   }
   void SetRef(uint * loc, const uint b)
   { uintSrc[b] = loc; dblSrc[b] = NULL; }

   void Sum(void);

   void Write(const ulong step, const bool firstPrint)
   {
      first = firstPrint;
      if (enable)
	 DoWrite(step);
   }

 private:

   std::string GetFName(std::string const& base, std::string const& uniqueName);

   void Zero(void)
   {
      for (uint b = 0; b < tot; b++)
	 block[b] = 0.0;
      samples = 0;
   }

   void DoWrite(const ulong step);

   bool first;
   std::ofstream outF;
   std::string name, varName;
   uint ** uintSrc, tot;
   double ** dblSrc;
   double * block, scl;
   uint samples;
   bool enable;
};

struct BlockAverages : OutputableBase
{
   BlockAverages(OutputVars & v){ this->var = &v; }

   ~BlockAverages(void) { if ( blocks != NULL ) delete[] blocks; }

   //No additional init.
   virtual void Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output);

   virtual void Sample(const ulong step);

   virtual void DoOutput(const ulong step);

 private:

   void InitVals(config_setup::EventSettings const& event)
   {
      stepsPerOut = event.frequency;
      invSteps = 1.0/stepsPerOut;
      enableOut = event.enable;
   }

   void AllocBlocks(void);

   void InitWatchSingle(config_setup::TrackedVars const& tracked);

   void InitWatchMulti(config_setup::TrackedVars const& tracked);

   //Block vars
   BlockAverage * blocks;
   uint numKindBlocks, totalBlocks;

   //Intermediate vars.
   uint samplesWrites;

   //Constants
   double invSteps;
};

#endif /*BLOCK_OUTPUT_H*/
