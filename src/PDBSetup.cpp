/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "PDBSetup.h" //Corresponding header to this body

#include <stdlib.h> //for exit

#include <map>    //for function lookup table.
#include <string> // for to_string
#include <vector>

#include "ConfigSetup.h"      //For restart info
#include "FixedWidthReader.h" //For fixed width reader
#include "MoveConst.h"
#include "StrLib.h" //for string comparison wrapper

#if BOX_TOTAL == 1
const std::string PDBSetup::pdbAlias[] = {"system PDB coordinate file"};
#else
const std::string PDBSetup::pdbAlias[] = {"box 0 PDB coordinate file",
                                          "box 1 PDB coordinate file"};
#endif

namespace pdb_setup {
void Remarks::SetRestart(config_setup::RestartSettings const &r) {
  restart = r.enable;
  recalcTrajectory = r.recalcTrajectory;
  restartFromXSC = r.restartFromXSCFile;
  restartFromBinary = r.restartFromBinaryCoorFile;

  for (uint b = 0; b < BOX_TOTAL; b++) {
    if (recalcTrajectory)
      reached[b] = false;
    else
      reached[b] = true;
  }
}
void Remarks::Read(FixedWidthReader &pdb) {
  using namespace pdb_entry::remark::field;
  using namespace pdb_entry::cryst1::field;

  if (restart) {
    // check if GOMC is taged and read the max dis, rot, vol value
    std::string varName;
    pdb.Get(varName, name::POS)
        .Get(disp[currBox], dis::POS)
        .Get(rotate[currBox], rot::POS)
        .Get(vol[currBox], vol::POS);

    CheckGOMC(varName);
  }
  if (recalcTrajectory) {
    std::string varName;
    pdb.Get(varName, name::POS)
        .Get(frameNumber[currBox], frameNum::POS)
        .Get(step[currBox], stepsNum::POS);

    if (frameNumber[currBox] == targetFrame[currBox])
      reached[currBox] = true;

    CheckGOMC(varName);
  }
}

void Remarks::CheckGOMC(std::string const &varName) {
  using namespace pdb_entry::remark::field;
  if (!str::compare(varName, name::STR_GOMC)) {
    std::cerr << "ERROR: " << "GOMC file's identifying tag "
              << "\"REMARK     GOMC\" is missing" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Remarks::Clear() { frameSteps.clear(); }

void Cryst1::Read(FixedWidthReader &pdb) {
  XYZ temp;
  using namespace pdb_entry::cryst1::field;
  hasVolume[currBox] = true;
  pdb.Get(temp.x, x::POS)
      .Get(temp.y, y::POS)
      .Get(temp.z, z::POS)
      .Get(cellAngle[currBox][0], ang_alpha::POS)
      .Get(cellAngle[currBox][1], ang_beta::POS)
      .Get(cellAngle[currBox][2], ang_gamma::POS);
  axis.Set(currBox, temp);
}

void Atoms::SetRestart(config_setup::RestartSettings const &r) {
  restart = r.enable;
  recalcTrajectory = r.recalcTrajectory;
}

void Atoms::Assign(std::string const &resName, const char l_chain,
                   const double l_x, const double l_y, const double l_z,
                   const double l_beta, const double l_occ) {
  // box.push_back((bool)(restart?(uint)(l_occ):currBox));
  beta.push_back(l_beta);
  occ.push_back(l_occ);
  box.push_back(currBox);
  ++numAtomsInBox[currBox];
  resNames.push_back(resName);
  chainLetter.push_back(l_chain);

  // push the coordinates of atoms to x, y, and z
  x.push_back(l_x);
  y.push_back(l_y);
  z.push_back(l_z);

  count++;
}

void Atoms::Read(FixedWidthReader &file) {
  using namespace pdb_entry::atom;
  char l_chain;
  uint resNum;
  std::string resName, atomName;
  double l_x, l_y, l_z, l_occ, l_beta;
  file.Get(atomName, field::alias::POS)
      .Get(resName, field::res_name::POS)
      .Get(resNum, field::res_num::POS)
      .Get(l_chain, field::chain::POS)
      .Get(l_x, field::x::POS)
      .Get(l_y, field::y::POS)
      .Get(l_z, field::z::POS)
      .Get(l_occ, field::occupancy::POS)
      .Get(l_beta, field::beta::POS);
  if (recalcTrajectory && (uint)l_occ != currBox) {
    return;
  }
  Assign(resName, l_chain, l_x, l_y, l_z, l_beta, l_occ);
}

void Atoms::Clear() {
  chainLetter.clear();
  x.clear();
  y.clear();
  z.clear();
  beta.clear();
  occ.clear();
  box.clear();
  resNames.clear();
  count = 0;
  for (uint b = 0; b < BOX_TOTAL; b++) {
    numAtomsInBox[b] = 0;
  }
  for (uint b = 0; b < BOX_TOTAL + 1; b++) {
    boxAtomOffset[0] = 0;
  }
}

// This method of finding minimum assumes all the box 0 atoms
// will be contiguous in the coordinates array.  This isn't
// the case on checkpoint restarts.  Since we went out of the
// way to ensure even when box transfers occur, the atoms
// remain in the same original order in the original data structure.
// We therefore can't rely on the molecule lookup to get the start
// and end of the box for the restarted data structures.
// Hence the numberOfAtoms array.

void Atoms::GetMinMaxAtoms(const uint b) {
  int stRange, endRange;

  boxAtomOffset[b + 1] = boxAtomOffset[b] + numAtomsInBox[b];

  // To prevent segfault, but we still need to set atomOffset
  if (numAtomsInBox[b] == 0)
    return;

  stRange = boxAtomOffset[b];
  endRange = boxAtomOffset[b + 1];

  min[b].x = *std::min_element(std::next(x.begin(), stRange),
                               std::next(x.begin(), endRange));
  min[b].y = *std::min_element(std::next(y.begin(), stRange),
                               std::next(y.begin(), endRange));
  min[b].z = *std::min_element(std::next(z.begin(), stRange),
                               std::next(z.begin(), endRange));
  max[b].x = *std::max_element(std::next(x.begin(), stRange),
                               std::next(x.begin(), endRange));
  max[b].y = *std::max_element(std::next(y.begin(), stRange),
                               std::next(y.begin(), endRange));
  max[b].z = *std::max_element(std::next(z.begin(), stRange),
                               std::next(z.begin(), endRange));
}

} // end namespace pdb_setup

void PDBSetup::Init(config_setup::RestartSettings const &restart,
                    std::string const *const name, uint frameNum) {
  // Clear the vectors for both atoms and remarks in case Init was called
  // more than once
  atoms.Clear();
  remarks.Clear();

  std::map<std::string, FWReadableBase *>::const_iterator dataKind;
  remarks.SetRestart(restart);
  atoms.SetRestart(restart);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    std::string varName = "";
    remarks.SetBox(b);
    remarks.SetFrameNumber(b, frameNum);
    cryst.SetBox(b);
    atoms.SetBox(b);
    std::string alias;
    if (remarks.recalcTrajectory) {
      sstrm::Converter toStr;
      std::string numStr = "";
      toStr << frameNum;
      toStr >> numStr;
      alias = pdbAlias[b] + " frame " + numStr;
    } else {
      alias = pdbAlias[b];
    }
    pdb[b].SetData(name[b], alias);

    // Open PDB only once and stay there
    // instead of re-opening it for every frame
    // refer to issue #131
    if (frameNum == 1)
      pdb[b].open();

    while (pdb[b].Read(varName, pdb_entry::label::POS)) {
      // If end of frame, and this is the frame we wanted,
      // end read on this file
      if (remarks.reached[b] && str::compare(varName, pdb_entry::end::STR)) {
        break;
      }

      // Call reader function if remarks were reached,
      // or it is a remark
      dataKind = dataKinds.find(varName);
      if (dataKind != dataKinds.end() &&
          (remarks.reached[b] ||
           str::compare(dataKind->first, pdb_entry::label::REMARK))) {
        dataKind->second->Read(pdb[b]);
      }
    }
    // If the recalcTrajectory is true and reached was still false
    // it means we couldn't find a remark and hence have to exit with error
    if (!remarks.reached[b] && remarks.recalcTrajectory) {
      std::cerr << "Error: Recalculate Trajectory is active..." << std::endl
                << ".. and couldn't find remark in PDB file!" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout.width(40);
    std::cout << std::left << "Finished reading: ";
    std::cout << "\t" << name[b] << std::endl;
    atoms.GetMinMaxAtoms(b);
  }
}

std::vector<ulong> PDBSetup::GetFrameSteps(std::string const *const name) {
  std::map<std::string, FWReadableBase *>::const_iterator dataKind;
  remarks.SetBox(mv::BOX0);
  FixedWidthReader pdb(name[mv::BOX0], pdbAlias[mv::BOX0]);
  pdb.open();
  std::string varName;
  uint count = 0;
  while (pdb.Read(varName, pdb_entry::label::POS)) {
    if (varName == pdb_entry::label::REMARK) {
      dataKind = dataKinds.find(varName);
      dataKind->second->Read(pdb);
      remarks.frameSteps.push_back(remarks.step[mv::BOX0]);
      count++;
    }
  }
  return remarks.frameSteps;
}
