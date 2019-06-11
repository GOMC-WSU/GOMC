/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PDB_OUTPUT_H
#define PDB_OUTPUT_H

#include <vector> //for molecule string storage.
#include <string> //to store lines of finished data.

#include "BasicTypes.h" //For uint

#include "OutputAbstracts.h"
#include "Molecules.h"
#include "MoleculeKind.h"
#include "StaticVals.h"
#include "Coordinates.h"
#include "Writer.h"
#include "PDBSetup.h" //For atoms class

class System;
namespace config_setup
{
struct Output;
}
class MoveSettings;
class MoleculeLookup;

struct PDBOutput : OutputableBase {
public:
  PDBOutput(System & sys, StaticVals const& statV);

  ~PDBOutput()
  {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      outF[b].close();
    }
  }

  //PDB does not need to sample on every step, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);

  virtual void DoOutput(const ulong step);
private:
  std::string GetDefaultAtomStr();

  void InitPartVec(pdb_setup::Atoms const& atoms);

  void SetMolBoxVec(std::vector<uint> & mBox);

  void PrintCryst1(const uint b, Writer & out);

  void PrintAtoms(const uint b, std::vector<uint> & mBox);

  //NEW_RESTART_CODE
  void DoOutputRebuildRestart(const uint step);
  void PrintAtomsRebuildRestart(const uint b);
  void PrintCrystRest(const uint b, const uint step, Writer & out);
  void PrintRemark(const uint b, const uint step, Writer & out);
  //NEW_RESTART_CODE

  void FormatAtom(std::string & line, const uint p, const uint m,
                  const char chain, std::string const& atomAlias,
                  std::string const& resName);

  void InsertAtomInLine(std::string & line, XYZ const& coor,
                        std::string const& occ, std::string const& beta);

  void PrintEnd(const uint b, Writer & out)
  {
    out.file << "END" << std::endl;
  }

  real ConvAng(const real t)
  {
    return acos(t) * 180.0 / M_PI;
  }

  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions& boxDimRef;
  Molecules const& molRef;
  Coordinates & coordCurrRef;
  COM & comCurrRef;

  Writer outF[BOX_TOTAL];
  //NEW_RESTART_CODE
  Writer outRebuildRestart[BOX_TOTAL];
  std::string outRebuildRestartFName[BOX_TOTAL];
  bool enableRestOut;
  ulong stepsRestPerOut;
  ulong stepsCoordPerOut;
  //NEW_RESTART_CODE
  bool enableOutState;
  std::vector<std::string> pStr;
  int frameNumber[BOX_TOTAL];
};

#endif /*PDB_OUTPUT_H*/
