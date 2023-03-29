/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DCLinear.h"

#include <cassert>

#include "DCOnSphere.h"
#include "DCRotateCOM.h"
#include "DCSingle.h"

using namespace cbmc;

DCLinear::DCLinear(System &sys, const Forcefield &ff, const MoleculeKind &kind,
                   const Setup &set)
    : data(sys, ff, set) {
  mol_setup::MolMap::const_iterator it = set.mol.kindMap.find(kind.uniqueName);
  assert(it != set.mol.kindMap.end());
  const mol_setup::MolKind setupKind = it->second;
  uint size = kind.NumAtoms();
  atomSize = size;

  idExchange = new DCRotateCOM(&data, setupKind);
  // First atom of the molecule
  forward.push_back(new DCSingle(&data, 0));
  backward.push_back(new DCSingle(&data, size - 1));
  // second atom of the molecule
  if (atomSize > 1) {
    forward.push_back(new DCOnSphere(&data, setupKind, 1, 0));
    backward.push_back(new DCOnSphere(&data, setupKind, size - 2, size - 1));
  }
}

DCLinear::~DCLinear() {
  for (uint i = 0; i < forward.size(); ++i) {
    delete forward[i];
    delete backward[i];
  }
  delete idExchange;
}

void DCLinear::Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  std::vector<DCComponent *> &comps = data.prng.randInt(1) ? forward : backward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }

  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // perform Intra-Swap move within the same box
  if (atomSize < 2) {
    Build(oldMol, newMol, molIndex);
    return;
  } else {
    // we only have two atoms in molecule: atom 0, 1
    uint fix = data.prng.randInt(1);
    // If fix == 0, forward (build atom 1), else backward (build atom 0)
    std::vector<DCComponent *> &comps = fix ? backward : forward;

    // copy the coordinate of the fix atom
    newMol.AddAtom(fix, oldMol.AtomPosition(fix));
    oldMol.ConfirmOldAtom(fix);
    // build the second atom
    comps[1]->PrepareNew(newMol, molIndex);
    comps[1]->BuildNew(newMol, molIndex);
    comps[1]->PrepareOld(oldMol, molIndex);
    comps[1]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // No crank shaft move for molecule with two atoms.
  // Instead we perform Regrowth move within the same box
  Regrowth(oldMol, newMol, molIndex);
}

void DCLinear::BuildIDNew(TrialMol &newMol, uint molIndex) {
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCLinear::BuildIDOld(TrialMol &oldMol, uint molIndex) {
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCLinear::BuildOld(TrialMol &oldMol, uint molIndex) {
  std::vector<DCComponent *> &comps = data.prng.randInt(1) ? forward : backward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::BuildNew(TrialMol &newMol, uint molIndex) {
  std::vector<DCComponent *> &comps = data.prng.randInt(1) ? forward : backward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }
}

void DCLinear::BuildGrowOld(TrialMol &oldMol, uint molIndex) {
  // If backbone is atom 0, we use forward, otherwise backward
  std::vector<DCComponent *> &comps = oldMol.GetAtomBB(0) ? backward : forward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}

void DCLinear::BuildGrowNew(TrialMol &newMol, uint molIndex) {
  // If backbone is atom 0, we use forward, otherwise backward
  std::vector<DCComponent *> &comps = newMol.GetAtomBB(0) ? backward : forward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
  }
}

void DCLinear::BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol,
                              uint molIndex) {
  // Get the seedIndex
  int sIndex;
  if (newMol.HasCav()) {
    sIndex = newMol.GetGrowingAtomIndex();
  } else if (oldMol.HasCav()) {
    sIndex = oldMol.GetGrowingAtomIndex();
  } else {
    std::cout << "Error: Calling BuildGrowInCav, but there is no cavity"
              << " defined for newMol and oldMol.\n";
    exit(EXIT_FAILURE);
  }

  // If backbone is atom 0, we use forward, otherwise backward
  std::vector<DCComponent *> &comps = sIndex ? backward : forward;
  for (uint i = 0; i < comps.size(); ++i) {
    comps[i]->PrepareNew(newMol, molIndex);
    comps[i]->BuildNew(newMol, molIndex);
    comps[i]->PrepareOld(oldMol, molIndex);
    comps[i]->BuildOld(oldMol, molIndex);
  }
}