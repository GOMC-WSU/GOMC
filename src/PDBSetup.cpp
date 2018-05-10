/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.30
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <vector>
#include <map> //for function lookup table.

#include "StrLib.h" //for string comparison wrapper
#include "PDBSetup.h" //Corresponding header to this body
#include "FixedWidthReader.h" //For fixed width reader
#include "ConfigSetup.h" //For restart info

#include <stdlib.h> //for exit

#if BOX_TOTAL == 1
const std::string PDBSetup::pdbAlias[] = {"system PDB coordinate file"};
#else
const std::string PDBSetup::pdbAlias[] = {"box 0 PDB coordinate file",
                                          "box 1 PDB coordinate file"
                                         };
#endif

namespace pdb_setup
{
void Remarks::SetRestart(config_setup::RestartSettings const& r )
{
  restart = r.enable;
  reached = true;
}
void Remarks::Read(FixedWidthReader & pdb)
{
  using namespace pdb_entry::remark::field;
  using namespace pdb_entry::cryst1::field;

  if(!restart)
    return;

  //check if GOMC is taged and read the max dis, rot, vol value
  std::string varName;
  pdb.Get(varName, name::POS)
  .Get(disp[currBox], dis::POS)
  .Get(rotate[currBox], rot::POS)
  .Get(vol[currBox], vol::POS);

  CheckGOMC(varName);
}

void Remarks::CheckGOMC(std::string const& varName)
{
  using namespace pdb_entry::remark::field;
  if (!str::compare(varName, name::STR_GOMC)) {
    std::cerr << "ERROR: Restart failed, "
              << "GOMC file's identifying tag "
              << "\"REMARK     GOMC\" is missing"
              << std::endl;
    exit(1);
  }
}

void Cryst1::Read(FixedWidthReader & pdb)
{
  XYZ temp;
  using namespace pdb_entry::cryst1::field;
  hasVolume = true;
  pdb.Get(temp.x, x::POS)
  .Get(temp.y, y::POS)
  .Get(temp.z, z::POS)
  .Get(cellAngle[currBox][0], ang_alpha::POS)
  .Get(cellAngle[currBox][1], ang_beta::POS)
  .Get(cellAngle[currBox][2], ang_gamma::POS);
  axis.Set(currBox, temp);
}

void Atoms::SetRestart(config_setup::RestartSettings const& r )
{
  restart = r.enable;
}

void Atoms::Assign(std::string const& atomName,
                   std::string const& resName,
                   const uint resNum,
                   const char l_chain, const double l_x,
                   const double l_y, const double l_z,
                   const double l_occ,
                   const double l_beta)
{
  //box.push_back((bool)(restart?(uint)(l_occ):currBox));
  beta.push_back(l_beta);
  box.push_back(currBox);
  atomAliases.push_back(atomName);
  resNamesFull.push_back(resName);
  if (resNum != currRes || resName != currResname || firstResInFile) {
    molBeta.push_back(l_beta);
    startIdxRes.push_back(count);
    currRes = resNum;
    currResname = resName;
    resNames.push_back(resName);
    chainLetter.push_back(l_chain);
    //Check if this kind of residue has been found
    uint kIndex = std::find(resKindNames.begin(), resKindNames.end(),
                            resName) - resKindNames.begin();
    // if not push it to resKindNames -> new molecule found
    if(kIndex == resKindNames.size()) {
      resKindNames.push_back(resName);
    }
    // pushes the index of the residue to the resKinds
    resKinds.push_back(kIndex);
  }
  // push the coordinates of atoms to x, y, and z
  x.push_back(l_x);
  y.push_back(l_y);
  z.push_back(l_z);

  count++;
  firstResInFile = false;
}

void Atoms::Read(FixedWidthReader & file)
{
  using namespace pdb_entry::atom;
  char l_chain;
  uint resNum;
  std::string resName, atomName;
  double l_x, l_y, l_z, l_occ, l_beta;
  file.Get(atomName, field::alias::POS)
  .Get(resName, field::res_name::POS)
  .Get(resNum, field::res_num::POS)
  .Get(l_chain, field::chain::POS).Get(l_x, field::x::POS)
  .Get(l_y, field::y::POS).Get(l_z, field::z::POS)
  .Get(l_occ, field::occupancy::POS)
  .Get(l_beta, field::beta::POS);
  Assign(atomName, resName, resNum, l_chain, l_x, l_y, l_z,
         l_occ, l_beta);
}

} //end namespace pdb_setup

void PDBSetup::Init(config_setup::RestartSettings const& restart,
                    std::string const*const name)
{
  using namespace std;
  map<string, FWReadableBase *>::const_iterator dataKind;
  remarks.SetRestart(restart);
  atoms.SetRestart(restart);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    std::string varName = "";
    remarks.SetBox(b);
    cryst.SetBox(b);
    atoms.SetBox(b);
    FixedWidthReader pdb(name[b], pdbAlias[b]);
    pdb.open();
    while (pdb.Read(varName, pdb_entry::label::POS)) {
      //If end of frame, and this is the frame we wanted,
      //end read on this file
      if (remarks.reached && str::compare(varName, pdb_entry::end::STR)) {
        break;
      }
      //Call reader function if remarks were reached,
      // or it is a remark
      dataKind = dataKinds.find(varName);
      if (dataKind != dataKinds.end() &&
          (remarks.reached ||
           str::compare(dataKind->first, pdb_entry::label::REMARK))) {
        dataKind->second->Read(pdb);
      }
    }
    pdb.close();
  }
}
