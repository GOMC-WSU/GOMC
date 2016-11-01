#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H

#include "BasicTypes.h" //For uint
#include "OutputAbstracts.h"
#include "Molecules.h"
#include "MoleculeKind.h"
#include "StaticVals.h"
#include "PDBSetup.h"
#include "MoveConst.h"
#include "OutputVars.h"

class System;
namespace config_setup { struct Output; }
class SystemPotential;
class Energy;
class Virial;
class MoveSettings;
class MoleculeLookup;

/*
struct CurrentAverages
{

   ~CurrentAverages()
   {
      if (molNum != NULL)
      {
	 delete[] molNum;
      }
      if (molFract != NULL)
      {
	    delete[] molFract;
      }
      if (density != NULL)
      {
	    delete[] density;
      }
   }

   void Init(const uint stEquil, const uint mkTot);

   void Sample(const uint step, double const*const lScale,
	       uint const*const lTries, uint const*const lAccept,
	       Energy const*const lEnergy, Virial const*const lVirial,
	       double const*const lVolume, double const*const molNum,
	       double const*const molFract, double const*const density);

   void OutputMoveInfo() const;
   void OutputEnAndVirInfo() const;
   void OutputBoxInfo() const;
   void OutputThermoInfo() const;

   uint stepsUntilEquil, numMolKinds;
   double scale[mv::SCALEABLE], tries[mv::COUNT], accepted[mv::COUNT];
   double preTries[mv::COUNT], preAccepted[mv::COUNT];
   Energy energy[BOXES_WITH_VOLUME];
   Virial virial[BOXES_WITH_VOLUME];
   double volume[BOXES_WITH_VOLUME], * molNum,
      * molFract, * density;
   double pressure, dHvDirect, dHvThermo;

};
*/

struct ConsoleOutput : OutputableBase
{
 public:
   ConsoleOutput(OutputVars & v){ this->var = &v; }

   //Console Output does not need to sample, so does nothing.
   virtual void Sample(const ulong step) {}
   
   virtual void Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output)
   {
      enableOut = output.console.enable;
      stepsPerOut = output.console.frequency;
      DoOutput(0);
   }

   virtual void DoOutput(const ulong step);

 private:

   void PrintBox(const uint box, const ulong step) const;

   void PrintSysStat() const;

   void PrintBanner(std::string const& str) const; 

   void PrintMoveKind(bool & somethingPrinted, const uint m,
		      const uint b, const ulong step) const;

   void PrintMoveStat(const uint box, const ulong step) const;

   void PrintMolKind(const uint k, const uint kb) const;

   void PrintEnergy(Energy const& en, Virial const& vir,
                    const bool intraOnly = false) const;
};

#endif /*CONSOLE_OUTPUT_H*/
