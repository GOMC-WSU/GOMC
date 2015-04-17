/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PDB_OUTPUT_H
#define PDB_OUTPUT_H

#include <vector> //for molecule string storage.
#include <string> //to store lines of finished data.

#include "../lib/BasicTypes.h" //For uint

#include "OutputAbstracts.h"
#include "Molecules.h"
#include "MoleculeKind.h"
#include "StaticVals.h"
#include "Coordinates.h"
#include "Writer.h"
#include "PDBSetup.h" //For atoms class

class System;
namespace config_setup { class Output; }
class MoveSettings;
class MoleculeLookup;

struct PDBOutput : OutputableBase
{
 public:
   PDBOutput(System & sys, StaticVals const& statV);

   //PDB does not need to sample on every step, so does nothing.
   virtual void Sample(const ulong step) {}

   virtual void Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output);

   virtual void DoOutput(const ulong step);
 private:
   std::string GetDefaultAtomStr();

   void InitPartVec(pdb_setup::Atoms const& atoms);

   void SetMolBoxVec(std::vector<uint> & mBox);

   void PrintCryst1(const uint b);

   void PrintAtoms(const uint b, std::vector<uint> & mBox);

   void PrintEnd(const uint b) { outF[b].file << "END" << std::endl; }
   
   MoveSettings & moveSetRef;
   MoleculeLookup & molLookupRef;
   BoxDimensions& boxDimRef;
   Molecules const& molRef;
   Coordinates & coordCurrRef;

   Writer outF[BOX_TOTAL];
   std::vector<std::string> pStr;
};

#endif /*PDB_OUTPUT_H*/

