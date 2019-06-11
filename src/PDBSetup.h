/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PDB_SETUP_H
#define PDB_SETUP_H

#include <vector>
#include <map> //for function lookup table.

#include "InputAbstracts.h" //For FWReadableBase
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "PDBConst.h" //For fields positions, etc.
#include "XYZArray.h" //For box dimensions.

namespace config_setup
{
struct RestartSettings;
}
struct FWReadableBase;


namespace pdb_setup
{
struct Remarks : FWReadableBase {
  uint currBox;
  real  disp[BOX_TOTAL], rotate[BOX_TOTAL], vol[BOX_TOTAL];
  ulong step[BOX_TOTAL];
  uint frameNumber[BOX_TOTAL], targetFrame[BOX_TOTAL];
  std::vector<ulong> frameSteps;
  bool restart, reached[BOX_TOTAL], recalcTrajectory;
  void SetRestart(config_setup::RestartSettings const& r);
  void Read(FixedWidthReader & pdb);
  void SetBox(const uint b)
  {
    currBox = b;
  }
  void SetFrameNumber(const uint b, const uint frameNum)
  {
    targetFrame[b] = frameNum;
  }
  void Clear();

private:
  void CheckGOMC(std::string const& varName);
};

struct Cryst1 : FWReadableBase {
  //box dimensions
  uint currBox;
  bool hasVolume;
  XYZArray axis;
  real cellAngle[BOX_TOTAL][3];
  Cryst1(void) : currBox(0), hasVolume(false), axis(BOX_TOTAL) {}
  void SetBox(const uint b)
  {
    currBox = b;
  }
  void Read(FixedWidthReader & pdb);
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
  }
  void Assign(std::string const& atomName,
              std::string const& resName,
              const uint resNum,
              const char l_chain,
              const real l_x,
              const real l_y,
              const real l_z,
              const real l_occ,
              const real l_beta);

  void Read(FixedWidthReader & file);
  void Clear();

  //private:
  //member data
  std::vector<char> chainLetter; //chain ids of each molecule
  std::vector<real> x, y, z; //coordinates of each particle
  std::vector<real> beta;  //beta value of each molecule
  std::vector<uint> box;
  std::vector<std::string> atomAliases, resNamesFull, resNames,
      resKindNames;
  std::vector<uint> startIdxRes, resKinds, molBeta;
  bool restart, firstResInFile, recalcTrajectory;
  //CurrRes is used to store res vals, currBox is used to
  //determine box either via the file (new) or the occupancy
  //(restart), count allows overwriting of coordinates during
  //second box read (restart only)
  uint currBox, count, currRes;
  std::string currResname;
};

}

struct PDBSetup {
  pdb_setup::Atoms atoms;
  pdb_setup::Cryst1 cryst;
  pdb_setup::Remarks remarks;
  PDBSetup(void) : dataKinds(SetReadFunctions()) {}
  void Init(config_setup::RestartSettings const& restart,
            std::string const*const name, uint frameNumber = 1);
  std::vector<ulong> GetFrameSteps(std::string const*const name);
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
