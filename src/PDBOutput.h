/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef PDB_OUTPUT_H
#define PDB_OUTPUT_H

#include <string> //to store lines of finished data.
#include <vector> //for molecule string storage.

#include "BasicTypes.h" //For uint
#include "Coordinates.h"
#include "MoleculeKind.h"
#include "Molecules.h"
#include "OutputAbstracts.h"
#include "PDBSetup.h" //For atoms class
#include "StaticVals.h"
#include "Writer.h"

class System;
namespace config_setup {
struct Output;
}
class MoveSettings;
class MoleculeLookup;

struct PDBOutput : OutputableBase {
public:
  PDBOutput(System &sys, StaticVals const &statV);

  ~PDBOutput() {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      outF[b].close();
    }
  }

  // PDB does not need to sample on every step, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  std::string GetDefaultAtomStr();

  void InitPartVec();

  void SetMolBoxVec(std::vector<uint> &mBox);

  void PrintCryst1(const uint b, Writer &out);

  void PrintAtoms(const uint b, std::vector<uint> &mBox);

  // NEW_RESTART_CODE
  void DoOutputRebuildRestart(const ulong step);
  void PrintAtomsRebuildRestart(const uint b);
  void PrintCrystRest(const uint b, const ulong step, Writer &out);
  void PrintRemark(const uint b, const ulong step, Writer &out);
  // NEW_RESTART_CODE

  void FormatAtom(std::string &line, const uint p, const uint m,
                  const char chain, std::string const &atomAlias,
                  std::string const &resName);

  template <typename T>
  void InsertAtomInLine(std::string &line, XYZ const &coor, const T &occ,
                        double const &beta);

  void PrintEnd(Writer &out) { out.file << "END" << std::endl; }

  double ConvAng(const double t) {
    // M_1_PI is 1/PI
    return acos(t) * 180.0 * M_1_PI;
  }

  MoveSettings &moveSetRef;
  MoleculeLookup &molLookupRef;
  BoxDimensions &boxDimRef;
  Molecules const &molRef;
  Coordinates &coordCurrRef;
  COM &comCurrRef;

  Writer outF[BOX_TOTAL];
  Writer outRebuildRestart[BOX_TOTAL];
  std::vector<std::string> pStr;
  int frameNumber[BOX_TOTAL];
};

#endif /*PDB_OUTPUT_H*/
