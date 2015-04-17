/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "PDBOutput.h"              //For spec;
#include "EnsemblePreprocessor.h"   //For BOX_TOTAL, ensemble
#include "System.h"                 //for init
#include "StaticVals.h"             //for init  
#include "MoleculeLookup.h"  //for lookup array (to get kind cnts, etc.)
#include "MoleculeKind.h"           //For kind names
#include "MoveSettings.h"           //For move settings/state
#include "PDBConst.h"               //For field locations/lengths

#include "../lib/StrStrmLib.h"      //For conversion from uint to string

#include <iostream>                 // for cout;

PDBOutput::PDBOutput(System & sys, StaticVals const& statV) :
   moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef), 
   coordCurrRef(sys.coordinates), 
   pStr(coordCurrRef.Count(),GetDefaultAtomStr()),
   boxDimRef(sys.boxDimRef), molRef(statV.mol) { }

std::string PDBOutput::GetDefaultAtomStr()
{
   using namespace pdb_entry::atom::field;
   using namespace pdb_entry;
   sstrm::Converter toStr;
   std::string defaultAtomStr(LINE_WIDTH, ' ');
   defaultAtomStr.replace(label::POS.START, label::POS.LENGTH, label::ATOM);
   return defaultAtomStr;
}

void PDBOutput::Init(pdb_setup::Atoms const& atoms,
		     config_setup::Output const& output)
{
   std::string bStr = "", aliasStr = "", numStr = "";
   sstrm::Converter toStr;
   enableOut = output.state.settings.enable;
   if (enableOut)
   {
      for (uint b = 0; b < BOX_TOTAL; ++b)
      {
	 //Get alias string, based on box #.
	 bStr = "Box ";
	 numStr = "";
	 toStr << b;
	 toStr >> numStr;
	 aliasStr = "Output PDB file for Box ";
	 aliasStr += bStr;
         bool notify;
#ifndef NDEBUG
         notify = true;
#else
         notify = false;
#endif
	 
         outF[b].Init(output.state.files.pdb.name[b], aliasStr, true, notify);
      }
      stepsPerOut = output.state.settings.frequency;
      InitPartVec(atoms);
      DoOutput(0);
   }
}

void PDBOutput::InitPartVec(pdb_setup::Atoms const& atoms)
{
   using namespace pdb_entry::atom::field;
   sstrm::Converter toStr;
   uint pStart = 0, pEnd = 0;
   //Start particle numbering @ 1
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b), 
	 end = molLookupRef.BoxEnd(b);
      while (m != end)
      {
	 uint mIndex = *m;
         std::string resName = atoms.resNames[mIndex];
	 char * chain = molRef.chain + mIndex;
         molRef.GetRangeStartStop(pStart, pEnd, mIndex);

         for (uint p = pStart; p < pEnd; ++p)
         {
	    std::string atomAlias = atoms.atomAliases[p];
            //Atom #
            toStr.Align(res_num::ALIGN).Replace(pStr[p], p + 1, atom_num::POS);

            uint posAliasStart = alias::POS.START;
            if (atomAlias.length() == 1)
            {
               ++posAliasStart;
            }
            pStr[p].replace(posAliasStart, alias::POS.LENGTH, atomAlias);

            //Res (molecule) name
            pStr[p].replace(res_name::POS.START, res_name::POS.LENGTH,
                            resName);
            //Res (molecule) chain (letter)
            pStr[p][chain::POS.START] = *chain;
            //Res (molecule) # -- add 1 to start counting @ 1
            toStr.Align(res_num::ALIGN).Replace(pStr[p], mIndex + 1,
                                                res_num::POS);
            toStr.Fixed().Align(beta::ALIGN).Precision(beta::PRECISION);
            toStr.Replace(pStr[p], beta::DEFAULT, beta::POS);
	 }
	 ++m;
      }
   }
}

void PDBOutput::DoOutput(const ulong step)
{     
   std::vector<uint> mBox(molRef.count);
   SetMolBoxVec(mBox);
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      outF[b].open();
      PrintCryst1(b);
      PrintAtoms(b, mBox);
      PrintEnd(b);
      outF[b].close();
   }
}

void PDBOutput::SetMolBoxVec(std::vector<uint> & mBox)
{
   for (uint b = 0; b < BOX_TOTAL; ++b)
   {
      MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b), 
	 end = molLookupRef.BoxEnd(b);
      while (m != end)
      {
	 mBox[*m] = b;
	 ++m;
      }
   }
}

void PDBOutput::PrintCryst1(const uint b)
{
   using namespace pdb_entry::cryst1::field;
   using namespace pdb_entry;
   sstrm::Converter toStr;
   std::string outStr(pdb_entry::LINE_WIDTH, ' '); 
   XYZ axis = boxDimRef.axis.Get(b);
   //Tag for crystallography -- cell dimensions.
   outStr.replace(label::POS.START, label::POS.LENGTH, label::CRYST1);
   //Add box dimensions
   toStr.Fixed().Align(x::ALIGN).Precision(x::PRECISION);
   toStr.Replace(outStr, axis.x, x::POS);
   toStr.Fixed().Align(y::ALIGN).Precision(y::PRECISION);
   toStr.Replace(outStr, axis.y, y::POS);
   toStr.Fixed().Align(z::ALIGN).Precision(z::PRECISION);
   toStr.Replace(outStr, axis.z, z::POS);
   //Add facet angles.
   toStr.Fixed().Align(ang_alpha::ALIGN).Precision(ang_alpha::PRECISION);
   toStr.Replace(outStr, ang_alpha::DEFAULT, ang_alpha::POS);
   toStr.Fixed().Align(ang_beta::ALIGN).Precision(ang_beta::PRECISION);
   toStr.Replace(outStr, ang_beta::DEFAULT, ang_beta::POS);
   toStr.Fixed().Align(ang_gamma::ALIGN).Precision(ang_gamma::PRECISION);
   toStr.Replace(outStr, ang_gamma::DEFAULT, ang_gamma::POS);
   //Add extra text junk.
   outStr.replace(space::POS.START, space::POS.LENGTH, space::DEFAULT);
   outStr.replace(zvalue::POS.START, zvalue::POS.LENGTH, zvalue::DEFAULT);
   //Write cell line
   outF[b].file << outStr << std::endl;
}

void PDBOutput::PrintAtoms(const uint b, std::vector<uint> & mBox)
{
   using namespace pdb_entry::atom::field;
   using namespace pdb_entry;
   sstrm::Converter toStr;
   bool inThisBox = false;
   uint pStart = 0, pEnd = 0;
   //Loop through all molecules
   for (uint m = 0; m < molRef.count; ++m)
   {
      //Loop through particles in mol.
      molRef.GetRangeStartStop(pStart, pEnd, m);
      inThisBox = (mBox[m]==b);
      for (uint p = pStart; p < pEnd; ++p)
      {
         XYZ coor;
         if (inThisBox)
         {
            coor = coordCurrRef.Get(p);
         }

	 //Fill in particle's stock string with new x, y, z, and occupancy
	 toStr.Fixed().Align(x::ALIGN).Precision(x::PRECISION);
	 toStr.Replace(pStr[p], coor.x, x::POS);
	 toStr.Fixed().Align(y::ALIGN).Precision(y::PRECISION);
	 toStr.Replace(pStr[p], coor.y, y::POS);
	 toStr.Fixed().Align(z::ALIGN).Precision(z::PRECISION);
	 toStr.Replace(pStr[p], coor.z, z::POS);
	 toStr.Align(occupancy::ALIGN);
	 toStr.Replace(pStr[p], occupancy::BOX[mBox[m]], occupancy::POS);
	 //Write finished string out.
	 outF[b].file << pStr[p] << std::endl;
      }
   }
}

