/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLECULES_H
#define MOLECULES_H

#include <map>
#include <string>

#include "BasicTypes.h" //For uint
#include "MolSetup.h"
// For iota
#include <numeric>
// For custom sort
#include "AlphaNum.h"

namespace pdb_setup {
class Atoms;
}
class FFSetup;
class Forcefield;
class System;

#include "MoleculeKind.h" //For member var.

// Note: This info is static and will never change in current ensembles
class Molecules {
public:
  Molecules();
  ~Molecules();
  bool operator==(const Molecules &other);

  const MoleculeKind &GetKind(const uint molIndex) const {
    return kinds[kIndex[molIndex]];
  }

  uint GetMolKind(const uint molIndex) const { return kIndex[molIndex]; }

  void Init(Setup &setup, Forcefield &forcefield, System &sys);

  uint NumAtomsByMol(const uint m) const { return start[m + 1] - start[m]; }
  uint NumAtoms(const uint mk) const { return kinds[mk].NumAtoms(); }

  int MolStart(const uint molIndex) const { return start[molIndex]; }

  int MolEnd(const uint molIndex) const { return start[molIndex + 1]; }

  int MolLength(const uint molIndex) const {
    return MolEnd(molIndex) - MolStart(molIndex);
  }

  void GetRange(uint &_start, uint &stop, uint &len, const uint m) const {
    _start = start[m];
    stop = start[m + 1];
    len = stop - _start;
  }

  void GetRangeStartStop(uint &_start, uint &stop, const uint m) const {
    _start = start[m];
    stop = start[m + 1];
  }

  void GetRestartOrderedRangeStartStop(uint &_start, uint &stop,
                                       const uint m) const {
    _start = restartOrderedStart[m];
    stop = restartOrderedStart[m + 1];
  }

  void GetRangeStartLength(uint &_start, uint &len, const uint m) const {
    _start = start[m];
    len = start[m + 1] - _start;
  }

  uint GetKindsCount() const { return kindsCount; }

  void PrintLJInfo(std::vector<uint> &totAtomKind,
                   std::vector<std::string> &names, Forcefield &forcefield);

  // private:
  // Kind index of each molecule and start in master particle array
  // Plus counts
  uint32_t *start;
  /* From checkpoint for loading binary coord/vel into the original PDBAtoms
   * object */
  uint32_t *restartOrderedStart;

  /* only used for output */
  uint32_t count;
  uint32_t atomCount;
  uint32_t *kIndex;
  uint32_t kIndexCount;
  uint *countByKind;
  char *chain;
  double *beta;
  double *occ;

  MoleculeKind *kinds;
  uint kindsCount;
  // These are never initialized or used
  // uint fractionKind;
  // uint lambdaSize;
  double *pairEnCorrections;
  double *pairVirCorrections;

  bool printFlag, restartFromCheckpoint;
};

#endif /*MOLECULES_H*/
