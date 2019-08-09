/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PSFOUTPUT_H
#define PSFOUTPUT_H

#include <vector>
#include <string>

#include "BasicTypes.h"
#include "MolSetup.h"
#include "System.h"
#include "MoleculeLookup.h"

class Molecules;

class PSFOutput
{
public:
  PSFOutput(const Molecules& molecules, const System &sys,
            mol_setup::MolMap& molMap,
            const std::vector<std::string>& kindNames);

  //Output PSF file to filename using default remarks
  void PrintPSF(const std::string& filename) const;

  //Output PSF file to filename, specifying remarks
  void PrintPSF(const std::string& filename, const std::vector<std::string>& remarks) const;

private:
  const Molecules* molecules;
  const MoleculeLookup & molLookRef;
  std::vector<mol_setup::MolKind> molKinds;
  std::vector<std::string> molNames;
  uint totalAtoms;
  uint totalBonds;
  uint totalAngles;
  uint totalDihs;

  void PrintRemarks(FILE* outfile, const std::vector<std::string>& remarks) const;
  void PrintAtoms(FILE* outfile) const;
  void PrintBonds(FILE* outfile) const;
  void PrintAngles(FILE* outfile) const;
  void PrintDihedrals(FILE* outfile) const;

  void CountMolecules();
};

#endif
