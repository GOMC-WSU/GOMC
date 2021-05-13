/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "PDBOutput.h"              //For spec;
#include "EnsemblePreprocessor.h"   //For BOX_TOTAL, ensemble
#include "System.h"                 //for init
#include "StaticVals.h"             //for init
#include "MoleculeLookup.h"         //for lookup array (to get kind cnts, etc.)
#include "MoleculeKind.h"           //For kind names
#include "MoveSettings.h"           //For move settings/state
#include "PDBConst.h"               //For field locations/lengths
#include "StrStrmLib.h"             //For conversion from uint to string
#include <iostream>                 //for cout;

PDBOutput::PDBOutput(System  & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  boxDimRef(sys.boxDimRef), molRef(statV.mol), coordCurrRef(sys.coordinates),
  comCurrRef(sys.com), pStr(coordCurrRef.Count(), GetDefaultAtomStr())
{
  for(int i = 0; i < BOX_TOTAL; i++)
    frameNumber[i] = 0;
}

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
  /* My reasoning is if you pass a PDB trajectory in, you want a PDB traj out 
  If we want to convert PDBTraj to DCDTraj use this version
  enableOut = output.state.settings.enable || atoms.recalcTrajectory;
  */
   /* For testing purposes it is good to print PDB since it is readable.
  enableOut = output.state.settings.enable || (atoms.recalcTrajectory && !atoms.recalcTrajectoryBinary);
  */
  enableOut = output.state.settings.enable || atoms.recalcTrajectory;
  enableRestOut = output.restart.settings.enable;

  stepsPerOut = output.state.settings.frequency;
  stepsRestPerOut = output.restart.settings.frequency;

  if (stepsPerOut > stepsRestPerOut)
    stepsPerOut = stepsRestPerOut;

  if (enableOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      //Get alias string, based on box #.
      bStr = "Box ";
      numStr = "";
      toStr << b + 1;
      toStr >> numStr;
      aliasStr = "Output PDB file for Box ";
      aliasStr += numStr;
      bool notify;
#ifndef NDEBUG
      notify = true;
#else
      notify = false;
#endif

      outF[b].Init(output.state.files.pdb.name[b], aliasStr, true, notify);
      outF[b].open();
    }
    InitPartVec();
    DoOutput(0);
  }

  if (enableRestOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      //Get alias string, based on box #.
      bStr = "Box ";
      numStr = "";
      toStr << b + 1;
      toStr >> numStr;
      aliasStr = "Output restart PDB file for Box ";
      aliasStr += numStr;
      bool notify;
#ifndef NDEBUG
      notify = true;
#else
      notify = false;
#endif
      outRebuildRestart[b].Init(output.restart.files.pdb.name[b], aliasStr,
                                true, notify);
    }
  }
}

void PDBOutput::InitPartVec()
{
  uint placementStart = 0, placementEnd = 0, dataStart = 0, dataEnd = 0, molecule = 0, atomIndex = 0, mI = 0;
  //Start particle numbering @ 1
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                 end = molLookupRef.BoxEnd(b);

    while (m != end) {
      if(molLookupRef.restartFromCheckpoint){
        mI = molLookupRef.constantReverseSort[molLookupRef.originalMoleculeIndices[molLookupRef.constantReverseSort[*m]]];
      } else {
        mI = *m;
      } 

      /* If we are on an original run, placement and data are equal.
      If we a checkpoint restarting they may differ
      */
      molRef.GetOriginalRangeStartStop(placementStart, placementEnd, mI);
      molRef.GetRangeStartStop(dataStart, dataEnd, *m);

      for (uint p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
        // If you don't want to preserve resID's comment this out -> mol = mI
        molecule = mI;
        if (molRef.kinds[molRef.originalKIndex[mI]].isMultiResidue){
          FormatAtom(pStr[p], p, molecule + molRef.kinds[molRef.originalKIndex[mI]].intraMoleculeResIDs[p - placementStart], 
                    molRef.chain[d],
                    molRef.kinds[molRef.originalKIndex[mI]].atomNames[p - placementStart], 
                    molRef.kinds[molRef.originalKIndex[mI]].resNames[p - placementStart]);
        } else {
          FormatAtom(pStr[p], p, molecule, 
                    molRef.chain[d],
                    molRef.kinds[molRef.originalKIndex[mI]].atomNames[p - placementStart], 
                    molRef.kinds[molRef.originalKIndex[mI]].resNames[p - placementStart]);
        }
        ++atomIndex;
      }
      ++m;
      ++molecule;
      /* If you want to keep orig resID's comment these out */
      if (molRef.kinds[molRef.originalKIndex[mI]].isMultiResidue){
        molecule += molRef.kinds[molRef.originalKIndex[mI]].intraMoleculeResIDs.back();
      }
      /* 0 & 9999 since FormatAtom adds 1 shifting to 1 and 10,000*/
      if(molecule == 9999)
        molecule = 0;
      /* If you want to keep orig resID's comment these out */
    }
  }
}

void PDBOutput::FormatAtom
(std::string & line, const uint p, const uint m, const char chain,
 std::string const& atomAlias, std::string const& resName)
{
  using namespace pdb_entry::atom::field;
  using namespace pdb_entry;
  sstrm::Converter toStr;
  //Atom #
  if(p + 1 < 100000) {
    toStr.Align(res_num::ALIGN).Replace(line, p + 1, atom_num::POS);
  } else {
    toStr.Align(res_num::ALIGN).Replace(line, "*****", atom_num::POS);
  }

  uint posAliasStart = alias::POS.START;
  if (atomAlias.length() == 1) {
    ++posAliasStart;
  }
  line.replace(posAliasStart, alias::POS.LENGTH, atomAlias);

  //Res (molecule) name
  line.replace(res_name::POS.START, res_name::POS.LENGTH,
               resName);
  //Res (molecule) chain (letter)
  line[chain::POS.START] = chain;
  //Res (molecule) # -- add 1 to start counting @ 1
  toStr.Align(res_num::ALIGN).Replace(line, (m % 9999) + 1, res_num::POS);

  toStr.Fixed().Align(beta::ALIGN).Precision(beta::PRECISION);
  toStr.Replace(line, beta::DEFAULT, beta::POS);
}

void PDBOutput::DoOutput(const ulong step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::PDB_OUTPUT);
  std::vector<uint> mBox(molRef.count);
  SetMolBoxVec(mBox);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    PrintRemark(b, step, outF[b]);
    PrintCryst1(b, outF[b]);
    PrintAtoms(b, mBox);
    PrintEnd(outF[b]);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PDB_OUTPUT);
}

//NEW_RESTART_CODE
void PDBOutput::DoOutputRestart(const ulong step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::PDB_RESTART_OUTPUT);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    outRebuildRestart[b].openOverwrite();
    PrintCrystRest(b, step, outRebuildRestart[b]);
    PrintCryst1(b, outRebuildRestart[b]);
    PrintAtomsRebuildRestart(b);
    PrintEnd(outRebuildRestart[b]);
    outRebuildRestart[b].close();
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PDB_RESTART_OUTPUT);
}

//NEW_RESTART_CODE
void PDBOutput::DoOutputRebuildRestart(const ulong step)
{
  GOMC_EVENT_START(1, GomcProfileEvent::PDB_RESTART_OUTPUT);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    outRebuildRestart[b].openOverwrite();
    PrintCrystRest(b, step, outRebuildRestart[b]);
    PrintCryst1(b, outRebuildRestart[b]);
    PrintAtomsRebuildRestart(b);
    PrintEnd(outRebuildRestart[b]);
    outRebuildRestart[b].close();
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PDB_RESTART_OUTPUT);
}
//NEW_RESTART_CODE

void PDBOutput::SetMolBoxVec(std::vector<uint> & mBox)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                 end = molLookupRef.BoxEnd(b);
    while (m != end) {
      mBox[*m] = b;
      ++m;
    }
  }
}

void PDBOutput::PrintCryst1(const uint b, Writer & out)
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
  toStr.Replace(outStr, ConvAng(boxDimRef.cosAngle[b][0]), ang_alpha::POS);
  toStr.Fixed().Align(ang_beta::ALIGN).Precision(ang_beta::PRECISION);
  toStr.Replace(outStr, ConvAng(boxDimRef.cosAngle[b][1]), ang_beta::POS);
  toStr.Fixed().Align(ang_gamma::ALIGN).Precision(ang_gamma::PRECISION);
  toStr.Replace(outStr, ConvAng(boxDimRef.cosAngle[b][2]), ang_gamma::POS);
  //Add extra text junk.
  outStr.replace(space::POS.START, space::POS.LENGTH, space::DEFAULT);
  outStr.replace(zvalue::POS.START, zvalue::POS.LENGTH, zvalue::DEFAULT);
  //Write cell line
  out.file << outStr << std::endl;
}

void PDBOutput::PrintCrystRest(const uint b, const ulong step, Writer & out)
{
  using namespace pdb_entry::cryst1::field;
  using namespace pdb_entry;
  using namespace pdb_entry::remark::field;
  double displace = moveSetRef.GetScaleTot(b, mv::DISPLACE);
  double rotate = moveSetRef.GetScaleTot(b, mv::ROTATE);
  double volume = 0.0;
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  volume = moveSetRef.GetScaleTot(b, mv::VOL_TRANSFER);
#endif
  sstrm::Converter toStr;
  std::string outStr(pdb_entry::LINE_WIDTH, ' ');
  //Tag for remark
  outStr.replace(label::POS.START, label::POS.LENGTH, label::REMARK);
  //Tag GOMC
  outStr.replace(name::POS.START, name::POS.LENGTH, name::STR_GOMC);
  //Add max amount of displacement, rotate, and volume
  toStr.Fixed().Align(dis::ALIGN).Precision(dis::PRECISION);
  toStr.Replace(outStr, displace, dis::POS);
  toStr.Fixed().Align(rot::ALIGN).Precision(rot::PRECISION);
  toStr.Replace(outStr, rotate, rot::POS);
  toStr.Fixed().Align(vol::ALIGN).Precision(vol::PRECISION);
  toStr.Replace(outStr, volume, vol::POS);
  //Add steps number
  toStr.Fixed().Align(stepsNum::ALIGN).Precision(stepsNum::PRECISION);
  toStr.Replace(outStr, step, stepsNum::POS);
  //Write cell line
  out.file << outStr << std::endl;
}


void PDBOutput::InsertAtomInLine(std::string & line, XYZ const& coor,
                                 std::string const& occ,
                                 std::string const& beta)
{
  using namespace pdb_entry::atom::field;
  using namespace pdb_entry;
  sstrm::Converter toStr;
  //Fill in particle's stock string with new x, y, z, and occupancy
  toStr.Fixed().Align(x::ALIGN).Precision(x::PRECISION);
  toStr.Replace(line, coor.x, x::POS);
  toStr.Fixed().Align(y::ALIGN).Precision(y::PRECISION);
  toStr.Replace(line, coor.y, y::POS);
  toStr.Fixed().Align(z::ALIGN).Precision(z::PRECISION);
  toStr.Replace(line, coor.z, z::POS);
  toStr.Align(occupancy::ALIGN);
  toStr.Replace(line, occ, occupancy::POS);
  toStr.Align(beta::ALIGN);
  toStr.Replace(line, beta, beta::POS);
}

void PDBOutput::PrintAtoms(const uint b, std::vector<uint> & mBox)
{
  using namespace pdb_entry::atom::field;
  using namespace pdb_entry;
  bool inThisBox = false;

  uint placementStart = 0, placementEnd = 0, dataStart = 0, dataEnd = 0, trajectoryI = 0, dataI = 0;
  //Start particle numbering @ 1
  for (uint box = 0; box < BOX_TOTAL; ++box) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(box),
                                  end = molLookupRef.BoxEnd(box);
    while (m != end) {
      if(molLookupRef.restartFromCheckpoint){
        trajectoryI = molLookupRef.constantReverseSort[molLookupRef.originalMoleculeIndices[molLookupRef.constantReverseSort[*m]]];
      } else {
        trajectoryI = *m;
      }       
      dataI = *m;
      //Loop through particles in mol.
      uint beta = molLookupRef.GetBeta(dataI);
      molRef.GetOriginalRangeStartStop(placementStart, placementEnd, trajectoryI);
      molRef.GetRangeStartStop(dataStart, dataEnd, dataI);
      XYZ ref = comCurrRef.Get(dataI);
      inThisBox = (mBox[dataI] == b);
      for (uint p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
        XYZ coor;
        if (inThisBox) {
          coor = coordCurrRef.Get(d);
          boxDimRef.UnwrapPBC(coor, b, ref);
        }
        InsertAtomInLine(pStr[p], coor, occupancy::BOX[mBox[dataI]], beta::FIX[beta]);
        //Write finished string out.
      }  
      ++m;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
    }
  }
  for (uint p = 0; p < coordCurrRef.Count(); ++p){
    outF[b].file << pStr[p] << std::endl;
  }
}

void PDBOutput::PrintAtomsRebuildRestart(const uint b)
{
  using namespace pdb_entry::atom::field;
  using namespace pdb_entry;
  uint molecule = 0, atom = 0, placementStart = 0, placementEnd = 0;
  for (uint k = 0; k < molRef.kindsCount; ++k) {
    uint countByKind = molLookupRef.NumKindInBox(k, b);
    for (uint kI = 0; kI < countByKind; ++kI) {
      uint molI = molLookupRef.GetMolNum(kI, k, b);
      /* We don't support hybrid fixed/flexible molecule currently */
      uint beta = molLookupRef.GetBeta(molI);
      molRef.GetRangeStartStop(placementStart, placementEnd, molI);
      XYZ ref = comCurrRef.Get(molI);
      for (uint p = placementStart; p < placementEnd; ++p) {
        std::string line = GetDefaultAtomStr();
        XYZ coor = coordCurrRef.Get(p);
        boxDimRef.UnwrapPBC(coor, b, ref);
        if (molRef.kinds[k].isMultiResidue){
          FormatAtom(line, atom, molecule + molRef.kinds[k].intraMoleculeResIDs[p - placementStart], molRef.chain[p],
                    molRef.kinds[k].atomNames[p - placementStart], molRef.kinds[k].resNames[p - placementStart]);
        } else {
          FormatAtom(line, atom, molecule, molRef.chain[p],
                    molRef.kinds[k].atomNames[p - placementStart], molRef.kinds[k].resNames[p - placementStart]);
        }
        //Fill in particle's stock string with new x, y, z, and occupancy
        InsertAtomInLine(line, coor, occupancy::BOX[0], beta::FIX[beta]);
        //Write finished string out.
        outRebuildRestart[b].file << line << std::endl;
        ++atom;
      }
      ++molecule;
      /* To add additional intramolecular residues */
      if (molRef.kinds[k].isMultiResidue){
        molecule += molRef.kinds[k].intraMoleculeResIDs.back();
      }
      /* 0 & 9999 since FormatAtom adds 1 shifting to 1 and 10,000*/
      if(molecule == 9999)
        molecule = 0;
    }
  }
}

// This function should print a remark for recalculating trajectory
// The format is the following:
// REMARK GOMC <frame number> <step number>
void PDBOutput::PrintRemark(const uint b, const ulong step, Writer & out)
{
  using namespace pdb_entry::cryst1::field;
  using namespace pdb_entry;
  using namespace pdb_entry::remark::field;

  sstrm::Converter toStr;
  std::string outStr(pdb_entry::LINE_WIDTH, ' ');
  //Tag for remark
  outStr.replace(label::POS.START, label::POS.LENGTH, label::REMARK);
  //Tag GOMC
  outStr.replace(name::POS.START, name::POS.LENGTH, name::STR_GOMC);

  // Print Frame number
  frameNumber[b]++;
  toStr.Fixed().Align(frameNum::ALIGN).Precision(frameNum::PRECISION);
  toStr.Replace(outStr, frameNumber[b], frameNum::POS);

  // Print step number
  toStr.Fixed().Align(stepsNum::ALIGN).Precision(stepsNum::PRECISION);
  toStr.Replace(outStr, step + 1, stepsNum::POS);

  //Write cell line
  out.file << outStr << std::endl;
}
