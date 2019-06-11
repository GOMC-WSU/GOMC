/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOLECULES_H
#define MOLECULES_H

#include "BasicTypes.h" //For uint
#include "MolSetup.h"
#include <map>
#include <string>

namespace pdb_setup
{
class Atoms;
}
class FFSetup;
class Forcefield;
class System;

#include "MoleculeKind.h" //For member var.

//Note: This info is static and will never change in current ensembles
class Molecules
{
public:
  Molecules();
  ~Molecules();

  const MoleculeKind& GetKind(const uint molIndex) const
  {
    return kinds[kIndex[molIndex]];
  }

  uint GetMolKind(const uint molIndex) const
  {
    return kIndex[molIndex];
  }

  void Init(Setup& setup, Forcefield& forcefield,
            System& sys);

  uint NumAtomsByMol(const uint m) const
  {
    return start[m + 1] - start[m];
  }
  uint NumAtoms(const uint mk) const
  {
    return kinds[mk].NumAtoms();
  }

  int MolStart(const uint molIndex) const
  {
    return start[molIndex];
  }

  int MolEnd(const uint molIndex) const
  {
    return start[molIndex + 1];
  }

  int MolLength(const uint molIndex) const
  {
    return MolEnd(molIndex) - MolStart(molIndex);
  }

  void GetRange(uint & _start, uint & stop, uint & len, const uint m) const
  {
    _start = start[m];
    stop = start[m + 1];
    len = stop - _start;
  }

  void GetRangeStartStop(uint & _start, uint & stop, const uint m) const
  {
    _start = start[m];
    stop = start[m + 1];
  }
  void GetRangeStartLength(uint & _start, uint & len, const uint m) const
  {
    _start = start[m];
    len = start[m + 1] - _start;
  }

  uint GetKindsCount() const
  {
    return kindsCount;
  }

  void PrintLJInfo(std::vector<uint> &totAtomKind,
                   std::vector<std::string> &names,
                   Forcefield & forcefield);

  //private:
  //Kind index of each molecule and start in master particle array
  //Plus counts
  uint* start;
  uint* kIndex;
  uint count;
  uint* countByKind;
  char* chain;

  MoleculeKind * kinds;
  uint kindsCount;
  real* pairEnCorrections;
  real* pairVirCorrections;

  bool printFlag;
};


#endif /*MOLECULES_H*/
