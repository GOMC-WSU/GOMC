/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "Geometry.h"

#include <algorithm>
#include <vector>

#include "FFSetup.h"
#include "MolSetup.h"

namespace {
// returns where p1 and p2 are found as a pair in bList, index+count if backward
int findPair(uint p1, uint p2, const BondList &bList) {
  for (uint i = 0; i < bList.count; ++i) {
    if (bList.part1[i] == p1 && bList.part2[i] == p2)
      return i;
    else if (bList.part1[i] == p2 && bList.part2[i] == p1)
      return i + bList.count;
  }
  return -1;
}
} // namespace

Nonbond::Nonbond() : part1(NULL), part2(NULL), count(0) {}

Nonbond::~Nonbond() {
  delete[] part1;
  delete[] part2;
}

void Nonbond::Init(const mol_setup::MolKind &molData) {
  unsigned int numAtoms = molData.atoms.size();
  // no nonbonded atoms
  if (numAtoms < 5) {
    return;
  }
  std::vector<unsigned int> part1Vec;
  std::vector<unsigned int> part2Vec;
  // because vector<bool> is weird
  std::vector<char> nonBondedAtoms(numAtoms);
  // find all the atoms that do not appear in an angle or dihedral
  // with the target atom (not checking bonds because either it's
  // diatomic or all bonds are part of angles
  // Also, ignore any index less than i to avoid double counting
  for (unsigned int i = 0; i < numAtoms; ++i) {
    nonBondedAtoms.assign(i, 0);
    nonBondedAtoms.insert(nonBondedAtoms.end(), numAtoms - i, 1);
    for (unsigned int j = 0; j < molData.angles.size(); ++j) {
      const mol_setup::Angle &ang = molData.angles[j];
      if (ang.a0 == i || ang.a1 == i || ang.a2 == i) {
        nonBondedAtoms[ang.a0] = 0;
        nonBondedAtoms[ang.a1] = 0;
        nonBondedAtoms[ang.a2] = 0;
      }
    }
    for (unsigned int j = 0; j < molData.dihedrals.size(); ++j) {
      const mol_setup::Dihedral &dih = molData.dihedrals[j];
      if (dih.a0 == i || dih.a1 == i || dih.a2 == i || dih.a3 == i) {
        nonBondedAtoms[dih.a0] = 0;
        nonBondedAtoms[dih.a1] = 0;
        nonBondedAtoms[dih.a2] = 0;
        nonBondedAtoms[dih.a3] = 0;
      }
    }
    // starting at i+1 to ignore double counting
    for (unsigned int j = i + 1; j < nonBondedAtoms.size(); ++j) {
      if (nonBondedAtoms[j] == 1) {
        part1Vec.push_back(i);
        part2Vec.push_back(j);
      }
    }
  }
  count = part1Vec.size();
  part1 = new uint[count];
  part2 = new uint[count];
  std::copy(part1Vec.begin(), part1Vec.end(), part1);
  std::copy(part2Vec.begin(), part2Vec.end(), part2);
}

void Nonbond_1_4::Init(const mol_setup::MolKind &molData) {
  unsigned int numAtoms = molData.atoms.size();
  // no nonbonded atoms
  if (numAtoms < 4) {
    return;
  }
  std::vector<unsigned int> part1Vec;
  std::vector<unsigned int> part2Vec;
  // because vector<bool> is weird
  std::vector<char> nonBondedAtoms(numAtoms);
  // find all the atoms that do not appear in an angle
  // with the target atom (not checking bonds because either it's
  // diatomic or all bonds are part of angles
  // Also, ignore any index less than i to avoid double counting
  for (unsigned int i = 0; i < numAtoms; ++i) {
    nonBondedAtoms.assign(i, 0);
    nonBondedAtoms.insert(nonBondedAtoms.end(), numAtoms - i, 0);
    for (unsigned int j = 0; j < molData.dihedrals.size(); ++j) {
      const mol_setup::Dihedral &dih = molData.dihedrals[j];
      if (dih.a0 == i || dih.a3 == i) {
        nonBondedAtoms[dih.a0] = 1;
        nonBondedAtoms[dih.a3] = 1;
      }
    }

    // starting at i+1 to ignore double counting
    for (unsigned int j = i + 1; j < nonBondedAtoms.size(); ++j) {
      if (nonBondedAtoms[j] == 1) {
        part1Vec.push_back(i);
        part2Vec.push_back(j);
      }
    }
  }
  count = part1Vec.size();
  part1 = new uint[count];
  part2 = new uint[count];
  std::copy(part1Vec.begin(), part1Vec.end(), part1);
  std::copy(part2Vec.begin(), part2Vec.end(), part2);
}

void Nonbond_1_3::Init(const mol_setup::MolKind &molData) {
  unsigned int numAtoms = molData.atoms.size();
  // no nonbonded atoms
  if (numAtoms < 3) {
    return;
  }
  std::vector<unsigned int> part1Vec;
  std::vector<unsigned int> part2Vec;
  // because vector<bool> is weird
  std::vector<char> nonBondedAtoms(numAtoms);
  // find all the atoms that form angle
  // with the target atom (not checking bonds because either it's
  // diatomic or all bonds are part of angles
  // Also, ignore any index less than i to avoid double counting
  for (unsigned int i = 0; i < numAtoms; ++i) {
    nonBondedAtoms.assign(i, 0);
    nonBondedAtoms.insert(nonBondedAtoms.end(), numAtoms - i, 0);

    for (unsigned int j = 0; j < molData.angles.size(); ++j) {
      const mol_setup::Angle &ang = molData.angles[j];
      if (ang.a0 == i || ang.a2 == i) {
        nonBondedAtoms[ang.a0] = 1;
        nonBondedAtoms[ang.a2] = 1;
      }
    }

    // starting at i+1 to ignore double counting
    for (unsigned int j = i + 1; j < nonBondedAtoms.size(); ++j) {
      if (nonBondedAtoms[j] == 1) {
        part1Vec.push_back(i);
        part2Vec.push_back(j);
      }
    }
  }
  count = part1Vec.size();
  part1 = new uint[count];
  part2 = new uint[count];
  std::copy(part1Vec.begin(), part1Vec.end(), part1);
  std::copy(part2Vec.begin(), part2Vec.end(), part2);
}

void EwaldNonbond::Init(const mol_setup::MolKind &molData) {
  unsigned int numAtoms = molData.atoms.size();
  std::vector<unsigned int> part1Vec;
  std::vector<unsigned int> part2Vec;
  std::vector<char> nonBondedAtoms(numAtoms);
  // find all possible pairs
  for (unsigned int i = 0; i < numAtoms; ++i) {
    // starting at i+1 to ignore double counting

    for (unsigned int j = i + 1; j < numAtoms; ++j) {
      part1Vec.push_back(i);
      part2Vec.push_back(j);
    }
  }
  count = part1Vec.size();
  part1 = new uint[count];
  part2 = new uint[count];
  std::copy(part1Vec.begin(), part1Vec.end(), part1);
  std::copy(part2Vec.begin(), part2Vec.end(), part2);
}

void BondList::Init(const std::vector<mol_setup::Bond> &bonds) {
  count = bonds.size();
  part1 = new uint[count];
  part2 = new uint[count];
  kinds = new uint[count];
  for (uint i = 0; i < count; ++i) {
    part1[i] = bonds[i].a0;
    part2[i] = bonds[i].a1;
    kinds[i] = bonds[i].kind;
  }
}

bool BondList::IsBonded(const uint &ai, const uint &aj) {
  for (uint i = 0; i < count; ++i) {
    if (part1[i] == ai && part2[i] == aj)
      return true;
    else if (part1[i] == aj && part2[i] == ai)
      return true;
  }
  return false;
}

BondList::BondList() : part1(NULL), part2(NULL), kinds(NULL), count(0) {}

BondList::~BondList() {
  delete[] part1;
  delete[] part2;
  delete[] kinds;
}

GeomFeature::GeomFeature(uint atomsPer)
    : bondsPer(atomsPer - 1), bondIndices(NULL), kinds(NULL), count(0) {}

GeomFeature::~GeomFeature() {
  delete[] bondIndices;
  delete[] kinds;
}

void GeomFeature::Init(const std::vector<mol_setup::Angle> &angles,
                       const BondList &bList) {
  count = angles.size();
  if (count == 0)
    return;
  // find corresponding bond indices
  kinds = new uint[count];
  bondIndices = new uint[count * bondsPer];

  int bondCounter = 0;
  for (uint i = 0; i < angles.size(); ++i) {
    bondIndices[bondCounter] = findPair(angles[i].a0, angles[i].a1, bList);
    ++bondCounter;
    bondIndices[bondCounter] = findPair(angles[i].a1, angles[i].a2, bList);
    ++bondCounter;
    kinds[i] = angles[i].kind;
  }
}

void GeomFeature::Init(const std::vector<mol_setup::Dihedral> &dihs,
                       const BondList &bList) {
  count = dihs.size();
  if (count == 0)
    return;
  // find corresponding bond indices
  kinds = new uint[count];
  bondIndices = new uint[count * bondsPer];

  int bondCounter = 0;
  for (uint i = 0; i < dihs.size(); ++i) {
    bondIndices[bondCounter] = findPair(dihs[i].a0, dihs[i].a1, bList);
    ++bondCounter;
    bondIndices[bondCounter] = findPair(dihs[i].a1, dihs[i].a2, bList);
    ++bondCounter;
    bondIndices[bondCounter] = findPair(dihs[i].a2, dihs[i].a3, bList);
    ++bondCounter;
    kinds[i] = dihs[i].kind;
  }
}

void GeomFeature::Init(const std::vector<mol_setup::Improper> &imps,
                       const BondList &bList) {
  count = imps.size();
  if (count == 0)
    return;
  // find corresponding bond indices
  kinds = new uint[count];
  bondIndices = new uint[count * bondsPer];

  int bondCounter = 0;
  for (uint i = 0; i < imps.size(); ++i) {
    bondIndices[bondCounter] = findPair(imps[i].a0, imps[i].a1, bList);
    ++bondCounter;
    bondIndices[bondCounter] = findPair(imps[i].a1, imps[i].a2, bList);
    ++bondCounter;
    bondIndices[bondCounter] = findPair(imps[i].a2, imps[i].a3, bList);
    ++bondCounter;
    kinds[i] = imps[i].kind;
  }
}

void SortedNonbond::Init(const Nonbond &nb, const uint numAtoms) {
  std::vector<uint> partnerVec;
  uint oldSize = 0;
  subdiv.Init(numAtoms);
  for (uint i = 0; i < numAtoms; ++i) {
    oldSize = partnerVec.size();
    for (uint j = 0; j < nb.count; ++j) {
      if (nb.part1[j] == i)
        partnerVec.push_back(nb.part2[j]);
      else if (nb.part2[j] == i)
        partnerVec.push_back(nb.part1[j]);
    }
    subdiv.Set(i, oldSize, partnerVec.size() - oldSize);
    // maybe sorting will improve access patterns?
    std::sort(partnerVec.begin() + oldSize, partnerVec.end());
  }
  partners = new uint[partnerVec.size()];
  std::copy(partnerVec.begin(), partnerVec.end(), partners);
}
