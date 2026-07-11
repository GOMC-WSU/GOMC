/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "MoleculeKind.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib> //for exit
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "CBMC.h"
#include "FFSetup.h"
#include "Forcefield.h"
#include "Geometry.h"
#include "MolSetup.h"
#include "PDBConst.h" //For resname length.
#include "PRNG.h"
#include "Setup.h"

void MoleculeKind::Init(uint &l_kindIndex, std::string const &l_name,
                        Setup const &setup, Forcefield const &forcefield,
                        System &sys) {
  kindIndex = l_kindIndex;
  mol_setup::MolMap::const_iterator dataIterator =
      setup.mol.kindMap.find(l_name);
  if (dataIterator == setup.mol.kindMap.end()) {
    std::cerr << "================================================" << std::endl
              << "Error: Molecule " << l_name
              << " not found in PDB file. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  const mol_setup::MolKind &molData = dataIterator->second;
  name = molData.moleculeName;
  uniqueName = l_name;
#if ENSEMBLE == GCMC
  std::map<std::string, double>::const_iterator
      kindCPIt = setup.config.sys.chemPot.cp.find(name),
      lastOne = setup.config.sys.chemPot.cp.end();

  // If we don't find a chemical potential for a kind in GCMC mode,
  // then quit.
  if (kindCPIt == lastOne) {
    std::cerr << "================================================" << std::endl
              << "Error: chemical potential is missing for "
              << molData.moleculeName << "." << std::endl
              << std::endl
              << "Here are the listed chemical potentials:" << std::endl
              << "----------------------------------------" << std::endl;

    // Print out whatever chemical potentials were read.
    for (kindCPIt = setup.config.sys.chemPot.cp.begin(); kindCPIt != lastOne;
         ++kindCPIt) {
      std::cerr << "Resname: " << kindCPIt->first
                << "      Value: " << kindCPIt->second << std::endl;
    }
    exit(EXIT_FAILURE);
  } else {
    chemPot = kindCPIt->second;
  }
#endif

  InitAtoms(molData);

  // Once-through topology objects

  oneThree = oneFour = false;

  if (forcefield.OneThree) {
    nonBonded_1_3.Init(molData);
    oneThree = true;
  }
  if (forcefield.OneFour) {
    nonBonded_1_4.Init(molData);
    oneFour = true;
  }
  if (forcefield.OneN) {
    nonBonded.Init(molData);
  }

  nonEwaldBonded.Init(molData);
  sortedNB_1_3.Init(nonBonded_1_3, numAtoms);
  sortedNB_1_4.Init(nonBonded_1_4, numAtoms);
  sortedNB.Init(nonBonded, numAtoms);
  sortedEwaldNB.Init(nonEwaldBonded, numAtoms);
  bondList.Init(molData.bonds);
  angles.Init(molData.angles, bondList);
  dihedrals.Init(molData.dihedrals, bondList);
  impropers.Init(molData.impropers, bondList);
  donorList.Init(molData.donors);
  acceptorList.Init(molData.acceptors);

#ifdef VARIABLE_PARTICLE_NUMBER
  builder = cbmc::MakeCBMC(sys, forcefield, *this, setup);
  // builder = new cbmc::LinearVlugt(sys, forcefield, *this, setup);
#endif
}

MoleculeKind::MoleculeKind()
    : angles(3), dihedrals(4), impropers(4), atomMass(NULL), builder(NULL),
      atomKind(NULL), atomCharge(NULL) {}

MoleculeKind::~MoleculeKind() {
  delete[] atomKind;
  delete[] atomMass;
  delete[] atomCharge;
  delete builder;
}

bool MoleculeKind::operator==(const MoleculeKind &other) {
  bool result = true;
  result &= (numAtoms == other.numAtoms);
  for (int i = 0; i < numAtoms; ++i) {
    result &= (atomKind[i] == other.atomKind[i]);
    result &= (atomMass[i] == other.atomMass[i]);
    result &= (atomCharge[i] == other.atomCharge[i]);
  }
  result &= (atomNames == other.atomNames);
  result &= (resNames == other.resNames);
  result &= (atomTypeNames == other.atomTypeNames);
  result &= (isMultiResidue == other.isMultiResidue);
  result &= (intraMoleculeResIDs == other.intraMoleculeResIDs);
  result &= (name == other.name);
  result &= (kindIndex == other.kindIndex);
  result &= (molMass == other.molMass);
#if ENSEMBLE == GCMC
  result &= (chemPot == other.chemPot);
#endif
  return result;
}

void MoleculeKind::InitAtoms(mol_setup::MolKind const &molData) {
  numAtoms = molData.atoms.size();
  atomKind = new uint[numAtoms];
  atomMass = new double[numAtoms];
  atomCharge = new double[numAtoms];
  molMass = 0;
  atomNames.clear();
  resNames.clear();

  /* These two entries all PSFOutput to
    correctly assign residueIDs to a map containing
    multi-residue and standard entries.  */
  isMultiResidue = molData.isMultiResidue;
  intraMoleculeResIDs = molData.intraMoleculeResIDs;

  // convert array of structures to structure of arrays
  for (uint i = 0; i < numAtoms; ++i) {
    const mol_setup::Atom &atom = molData.atoms[i];
    atomNames.push_back(atom.name);
    resNames.push_back(atom.residue);
    atomTypeNames.push_back(atom.type);
    atomMass[i] = atom.mass;
    molMass += atom.mass;
    atomCharge[i] = atom.charge;
    atomKind[i] = atom.kind;
  }
}

double MoleculeKind::GetMoleculeCharge() {
  double netCharge = 0.0;

  for (uint i = 0; i < numAtoms; ++i) {
    // to calculate net charge
    netCharge += atomCharge[i];
  }

  return netCharge;
}

bool MoleculeKind::MoleculeHasCharge() {
  for (uint i = 0; i < numAtoms; ++i) {
    if (std::abs(atomCharge[i]) > 0.0)
      return true;
  }
  return false;
}
