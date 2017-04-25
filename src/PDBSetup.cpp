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
const std::string PDBSetup::pdbAlias[] = {"Box 1 PDB coordinate file",
                                          "Box 2 PDB coordinate file"
                                         };
#endif

namespace pdb_setup
{
void Remarks::SetRestart(config_setup::RestartSettings const& r )
{
  restart = r.enable;
  reached = (!restart);
  restartStep = r.step;
}
void Remarks::Read(FixedWidthReader & pdb)
{
  using namespace pdb_entry::remark::field;
  if (!restart) return;
  uint remNum;
  ulong readStep;
  std::string varName;
  pdb.Get(remNum, rem_num::POS).Get(varName, name::POS)
  .Get(readStep, data::POS);
  HandleRemark(remNum, varName, readStep);
}

void Remarks::HandleRemark(const uint num,
                           std::string const& varName,
                           const ulong step)
{
  switch (num)
  {
  case 2:
    CheckStep(varName, step);
    break;
  case 1:
    CheckGOMC(varName);
  default:
    break;
  }
}

void Remarks::CheckStep(std::string const& varName,
                        const ulong readStep)
{
  using namespace pdb_entry::remark::field;
  if (!str::compare(varName, name::STR_STEP))   //malformed PDB
  {
    std::cerr << "ERROR: Restart failed, "
              << "GOMC file's step REMARK is "
              << "malformed." << std::endl;
    exit(1);
  }
  reached = (readStep == restartStep);
#ifndef NDEBUG
  if (reached && restart)
    std::cout << "Restart step " << restartStep << " reached."
              << std::endl;
#endif
}

void Remarks::CheckGOMC(std::string const& varName)
{
  using namespace pdb_entry::remark::field;
  if (!str::compare(varName, name::STR_GOMC))
  {
    std::cerr << "ERROR: Restart failed, "
              << "GOMC file's identifying tag "
              << "\"REMARK  1   GOMC\" is missing"
              << std::endl;
    exit(1);
  }
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
  if (!restart || currBox == 0)
  {
    //box.push_back((bool)(restart?(uint)(l_occ):currBox));
    beta.push_back(l_beta);
    box.push_back(currBox);
    atomAliases.push_back(atomName);
    resNamesFull.push_back(resName);
    if (resNum != currRes || firstResInFile)
    {
      startIdxRes.push_back(count);
      currRes = resNum;
      resNames.push_back(resName);
      chainLetter.push_back(l_chain);
      //Check if this kind of residue has been found
      uint kIndex = std::find(resKindNames.begin(),
                              resKindNames.end(),
                              resName) - resKindNames.begin();
      if (kIndex==resKindNames.size())
      {
        resKindNames.push_back(resName);
      }
      resKinds.push_back(kIndex);
    }
    x.push_back(l_x);
    y.push_back(l_y);
    z.push_back(l_z);
  }
  else if (box[count]==currBox)
  {
    //Overwrite members in 2nd box for restart file
    chainLetter[count] = l_chain;
    x[count] = l_x;
    y[count] = l_y;
    z[count] = l_z;
  }
  count++;
  firstResInFile = false;
}

void Atoms::Read(FixedWidthReader & file)
{
  using namespace pdb_entry::atom;
  char l_chain;
  uint resNum;
  std::string resName, atomName;
  double l_x, l_y, l_z, l_occ;
  file.Get(atomName, field::alias::POS)
  .Get(resName, field::res_name::POS)
  .Get(resNum, field::res_num::POS)
  .Get(l_chain,field::chain::POS).Get(l_x,field::x::POS)
  .Get(l_y,field::y::POS).Get(l_z,field::z::POS)
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
#ifndef NDEBUG
  std::cout << (restart.enable?
                "Reading GOMC dumped restart PDB file(s).":
                "Reading new system from PDB file(s).") << std::endl;
#endif
  for (uint b = 0; b < BOX_TOTAL; b++)
  {
    std::string varName="";
    cryst.SetBox(b);
    atoms.SetBox(b);
    FixedWidthReader pdb(name[b], pdbAlias[b]);
    pdb.open();
    while (pdb.Read(varName, pdb_entry::label::POS))
    {
      //If end of frame, and this is the frame we wanted,
      //end read on this
      //file
      if (remarks.reached &&
          str::compare(varName, pdb_entry::end::STR))
      {
        break;
      }
      //Call reader function if remarks were reached,
      // or it is a remark
      dataKind = dataKinds.find(varName);
      if (dataKind != dataKinds.end() &&
          (remarks.reached ||
           str::compare(dataKind->first,pdb_entry::label::REMARK)))
      {
        dataKind->second->Read(pdb);
      }
    }
    pdb.close();
  }
}
