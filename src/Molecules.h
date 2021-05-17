/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
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
// For iota
#include <numeric>
// For custom sort 
#include "AlphaNum.h"


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

  void GetOriginalRangeStartStop(uint & _start, uint & stop, const uint m) const
  {
    _start = originalStart[m];
    stop = originalStart[m + 1];
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
  uint32_t* start;
  /* From checkpoint for consistent trajectories */
  uint32_t* originalStart;
  uint32_t* originalKIndex;

  /* only used for output */
  uint32_t count;
  uint32_t* kIndex;
  uint32_t kIndexCount;
  uint* countByKind;
  char* chain;
  double* beta;


  MoleculeKind * kinds;
  uint kindsCount;
  uint fractionKind, lambdaSize;
  double* pairEnCorrections;
  double* pairVirCorrections;

  bool printFlag;
};


#endif /*MOLECULES_H*/
