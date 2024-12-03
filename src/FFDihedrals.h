/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FF_DIHEDRALS_H
#define FF_DIHEDRALS_H

#include "BasicTypes.h"      //For "uint"
#include "FFConst.h"         //GetRot
#include "NumLib.h"          //Sq
#include "SubdividedArray.h" //Subdivisions in master array for dih sets.

namespace ff_setup {
class Dihedral;
}

// FFDihKind
// Stores parameters for dihedral kinds
// dihedrals may require several sets of parameters to calculate,
// these are stored in consecutive array elements, each kind's
// parameter sets are delinated by subdiv
// Dihedrals may be calculated by several different styles
class FFDihedrals {
public:
  // calculate the energy of dih kind at angle phi
  double Calc(const uint kind, const double phi) const;
  // Initialize with data from parameter files
  void Init(ff_setup::Dihedral const &dih);

  FFDihedrals(void) : Kchi(NULL), delta(NULL), n(NULL) {}
  ~FFDihedrals(void);

private:
  // dih kind params
  SubdividedArray subdiv;
  double *Kchi, *delta;
  uint *n;
};

inline double FFDihedrals::Calc(const uint kind, const double phi) const {
  double sum = 0.0;
  for (uint i = subdiv.Begin(kind); i != subdiv.End(kind); ++i) {
    sum += Kchi[i] * (1 + cos(n[i] * phi - delta[i]));
  }
  return sum;
}

#endif /*FF_DIHEDRALS_H*/
