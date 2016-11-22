#ifndef PDB_SETUP_H
#define PDB_SETUP_H

#include <vector>
#include <map> //for function lookup table.

#include "InputAbstracts.h" //For FWReadableBase
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "PDBConst.h" //For fields positions, etc.
#include "XYZArray.h" //For box dimensions.

namespace config_setup { struct RestartSettings; }
struct FWReadableBase;


namespace pdb_setup
{
   struct Remarks : FWReadableBase
   {
      bool restart, reached;
      ulong restartStep;

      void SetRestart(config_setup::RestartSettings const& r);
      void Read(FixedWidthReader & pdb);

   private:
      void HandleRemark(const uint num, std::string const& varName,
         const ulong step);
      void CheckStep(std::string const& varName,
		     const ulong readStep);
      void CheckGOMC(std::string const& varName);
   };

   struct Cryst1 : FWReadableBase
   {
      //box dimensions
      uint currBox;
      bool hasVolume;
      XYZArray axis;
      Cryst1(void) : currBox(0), hasVolume(false), axis(BOX_TOTAL) {}
      void SetBox(const uint b) { currBox = b; }
      void Read(FixedWidthReader & pdb)
      {
         XYZ temp;
         using namespace pdb_entry::cryst1::field;
         hasVolume = true;
         pdb.Get(temp.x, x::POS)
	    .Get(temp.y, y::POS)
	    .Get(temp.z, z::POS);
         axis.Set(currBox, temp);
      }

   };

   class Atoms : public FWReadableBase
   {
   public:
      //Set the current residue to something other than 1
      Atoms(void) : restart(false), currBox(0), count(0),
	 currRes(10) {}
      void SetRestart(config_setup::RestartSettings const& r);
      void SetBox(const uint b)
      {
         currBox = b;
	 firstResInFile = true;
         //restart count if new system, second box. 
         count = ((b == 1 && restart) ? 0 : count);
      }
      void Assign(std::string const& atomName,
		  std::string const& resName,
		  const uint resNum,
		  const char l_chain,
		  const double l_x,
		  const double l_y,
		  const double l_z, 
		  const double l_occ);
      void Read(FixedWidthReader & file);

   //private:
	   //member data
	   std::vector<char> chainLetter; //chain ids of each molecule
	   std::vector<double> x, y, z; //coordinates of each particle
	   std::vector<uint> box;
	   std::vector<std::string> atomAliases, resNamesFull, resNames,
		   resKindNames;
	   std::vector<uint> startIdxRes, resKinds;
	   bool restart, firstResInFile;
	   //CurrRes is used to store res vals, currBox is used to 
	   //determine box either via the file (new) or the occupancy
	   //(restart), count allows overwriting of coordinates during
	   //second box read (restart only)
	   uint currBox, count, currRes;
   };

}

struct PDBSetup
{
   pdb_setup::Atoms atoms;
   pdb_setup::Cryst1 cryst;
   pdb_setup::Remarks remarks;
   PDBSetup(void) : dataKinds(SetReadFunctions()) {}
   void Init(config_setup::RestartSettings const& restart,
      std::string const*const name);
 private:
   //Map variable names to functions
   std::map<std::string, FWReadableBase *>  SetReadFunctions(void)
   {
      std::map<std::string, FWReadableBase *> funct;
      funct[pdb_entry::label::REMARK] = &remarks;
      funct[pdb_entry::label::CRYST1] = &cryst;
      funct[pdb_entry::label::ATOM] = &atoms;
      return funct;
   }
   const std::map<std::string, FWReadableBase *> dataKinds;
   static const std::string pdbAlias[];
};

#endif /*PDB_SETUP_H*/
