/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "CBMC.h"

#include <vector>

#include "DCCyclic.h"
#include "DCGraph.h"
#include "DCLinear.h"
#include "MolSetup.h"
#include "MoleculeKind.h"

namespace cbmc {

CBMC *MakeCBMC(System &sys, const Forcefield &ff, const MoleculeKind &kind,
               const Setup &set) {
  std::vector<uint> bondCount(kind.NumAtoms(), 0);
  for (uint i = 0; i < kind.bondList.count; ++i) {
    bondCount[kind.bondList.part1[i]]++;
    bondCount[kind.bondList.part2[i]]++;
  }

  bool cyclic = (kind.NumBonds() > kind.NumAtoms() - 1) ? true : false;

  if (cyclic) {
    return new DCCyclic(sys, ff, kind, set);
  } else if (kind.NumAtoms() > 2) {
    // Any molecule woth 3 atoms and more will be built in DCGraph
    return new DCGraph(sys, ff, kind, set);
  } else {
    return new DCLinear(sys, ff, kind, set);
  }
}

} // namespace cbmc
