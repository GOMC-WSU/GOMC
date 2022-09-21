/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLSETUP_H
#define MOLSETUP_H

#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <map>
#include <string>
#include <vector>

#include "AlphaNum.h"
#include "BasicTypes.h"
#include "BondAdjacencyList.h"
#include "EnsemblePreprocessor.h"

namespace config_setup {
struct RestartSettings;
}
namespace pdb_setup {
class Atoms;
}
class FFSetup;

namespace mol_setup {
struct MoleculeVariables {
  MoleculeVariables() = default;

  std::vector<uint> startIdxMolecules, moleculeKinds;
  // moleculeNames - length number of molecules
  // moleculeKindNames - length number of kinds
  // uniqueMapKeys - length number of kinds
  std::vector<std::string> moleculeNames, moleculeKindNames,
      moleculeSegmentNames, uniqueMapKeys;
  uint lastAtomIndexInBox0 = 0;
  uint numberMolsInBox0 = 0;
  uint molKindIndex = 0;
  uint stringSuffixMultiResidue = 0;
  uint stringSuffixNonMultiResidue = 0;
  uint moleculeIteration = 0;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &startIdxMolecules;
    ar &moleculeKinds;
    ar &moleculeNames;
    ar &moleculeKindNames;
    ar &moleculeSegmentNames;
    ar &uniqueMapKeys;
    ar &lastAtomIndexInBox0;
    ar &numberMolsInBox0;
    ar &molKindIndex;
    ar &stringSuffixMultiResidue;
    ar &stringSuffixNonMultiResidue;
    ar &moleculeIteration;
  }
};

//! structure to contain an atom's data during initialization
class Atom {
public:
  // Default constructor for serialization/deserialization
  Atom(void) = default;
  Atom(std::string const &l_name, std::string const &l_residue, uint l_resID,
       std::string const &l_segment, std::string const &l_type,
       const double l_charge, const double l_mass)
      : name(l_name), type(l_type), residue(l_residue), segment(l_segment),
        charge(l_charge), mass(l_mass), residueID(l_resID) {}
  // private:
  // name (within a molecule) and type (for forcefield params)
  std::string name, type, residue, segment;
  double charge, mass;
  // kind index
  /* ResID is by the PSF Parser to determine multi-residue status by comparing
    1st vs all of resIDs in a row in the moleculeXAtomIDY 2D vector */
  uint residueID;

  uint kind;

  bool operator==(const Atom &atm) const {
    if (type == atm.type && charge == atm.charge && residue == atm.residue &&
        mass == atm.mass)
      return true;
    else
      return false;
  }

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &name;
    ar &type;
    ar &residue;
    ar &segment;
    ar &charge;
    ar &mass;
    ar &residueID;
    ar &kind;
  }
};

class Dihedral {
public:
  // Default constructor for serialization/deserialization
  Dihedral(void) = default;
  Dihedral(uint atom0, uint atom1, uint atom2, uint atom3)
      : a0(atom0), a1(atom1), a2(atom2), a3(atom3) {}
  // some xplor PSF files have duplicate dihedrals, we need to ignore these
  bool operator==(const Dihedral &other) const;
  bool operator!=(const Dihedral &other) const;

  // private:
  // atoms
  uint a0, a1, a2, a3;
  uint kind;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &a0;
    ar &a1;
    ar &a2;
    ar &a3;
    ar &kind;
  }
};

class Improper {
public:
  // Default constructor for serialization/deserialization
  Improper(void) = default;
  Improper(uint atom0, uint atom1, uint atom2, uint atom3)
      : a0(atom0), a1(atom1), a2(atom2), a3(atom3) {}
  // Impropers are order specific, so these methods different from Dihedrals,
  // Hence why we define a separate class for Impropers
  bool operator==(const Improper &other) const;
  bool operator!=(const Improper &other) const;

  // private:
  // atoms
  uint a0, a1, a2, a3;
  uint kind;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &a0;
    ar &a1;
    ar &a2;
    ar &a3;
    ar &kind;
  }
};

class Angle {
public:
  // Default constructor for serialization/deserialization
  Angle(void) = default;
  Angle(uint atom0, uint atom1, uint atom2) : a0(atom0), a1(atom1), a2(atom2) {}

  // private:
  uint a0, a1, a2;
  uint kind;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &a0;
    ar &a1;
    ar &a2;
    ar &kind;
  }
};

class Bond {
public:
  // Default constructor for serialization/deserialization
  Bond(void) = default;
  Bond(uint atom0, uint atom1) : a0(atom0), a1(atom1) {}
  //   private:
  uint a0, a1;
  uint kind;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &a0;
    ar &a1;
    ar &kind;
  }
};

//! Structure to contain a molecule kind's data during initialization
class MolKind {
public:
  // Default constructor for serialization/deserialization
  MolKind() = default;

  // private:
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  std::vector<Improper> impropers;
  /*Technically H-Bonds, but no difference in class structure */
  std::vector<Bond> donors;
  std::vector<Bond> acceptors;

  uint kindIndex;

  // Used to search PSF file for geometry, meaningless after that
  uint firstAtomID, firstMolID;
  // true while the molecule is still open for modification during PSF read
  bool incomplete = true;
  bool isMultiResidue;
  std::vector<uint> intraMoleculeResIDs;
  // Used to map chemical potentials in config file to molecules
  // Note for proteins there is some guesswork.  As they are encountered
  // in the psf file they are named PROTA, PROTB, ..
  std::string moleculeName;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &atoms;
    ar &bonds;
    ar &angles;
    ar &dihedrals;
    ar &impropers;
    /*Technically H-Bonds, but no difference in class structure */
    ar &donors;
    ar &acceptors;

    ar &kindIndex;

    // Used to search PSF file for geometry, meaningless after that
    ar &firstAtomID;
    ar &firstMolID;
    // true while the molecule is still open for modification during PSF read
    ar &incomplete;
    ar &isMultiResidue;
    ar &intraMoleculeResIDs;
    ar &moleculeName;
  }
};

// List of dihedrals with atom at one end, atom first
std::vector<Dihedral> AtomEndDihs(const MolKind &molKind, uint atom);
// List of dihedrals with atom and partner in middle, atom in a1
std::vector<Dihedral> DihsOnBond(const MolKind &molKind, uint atom,
                                 uint partner);
// List of all dihedrals in the molecule kind
std::vector<Dihedral> DihsAll(const MolKind &molKind);
// List of angles with atom at one end, atom first
std::vector<Angle> AtomEndAngles(const MolKind &molKind, uint atom);
// List of angles with atom in middle
std::vector<Angle> AtomMidAngles(const MolKind &molKind, uint atom);
// List of angles with atom at one end, and mid in middle, atom first
std::vector<Angle> AtomMidEndAngles(const MolKind &molKind, uint mid,
                                    uint atom);
// List of all angles in the molecule kind
std::vector<Angle> AngsAll(const MolKind &molKind);
// List of bonds with atom at one end, atom first
std::vector<Bond> AtomBonds(const MolKind &molKind, uint atom);
// List of all bonds in the molecule kind
std::vector<Bond> BondsAll(const MolKind &molKind);

// first element (string) is name of molecule type
typedef std::map<std::string, MolKind> MolMap;
typedef std::map<std::size_t, std::vector<std::string>> SizeMap;

//! Reads one or more PSF files into kindMap
/*!
 *\param kindMap map to add PSF data to
 *\param psfFilename array of strings containing filenames
 *\param numFiles number of files to read
 *\return -1 if failed, 0 if successful
 */
int ReadCombinePSF(MoleculeVariables &molVars, MolMap &kindMap,
                   SizeMap &sizeMap, const std::string *psfFilename,
                   const bool *psfDefined, pdb_setup::Atoms &pdbAtoms);

void PrintMolMapVerbose(const MolMap &kindMap);
void PrintMolMapBrief(const MolMap &kindMap);
} // namespace mol_setup

// wrapper struct for consistent interface
class MolSetup {
public:
  MolSetup() = default;
  class Atom;
  void createKindMap(mol_setup::MoleculeVariables &molVars,
                     const BondAdjacencyList &bondAdjList,
                     const std::vector<std::vector<uint>> &moleculeXAtomIDY,
                     std::vector<mol_setup::Atom> &allAtoms,
                     mol_setup::MolMap &kindMap, mol_setup::SizeMap &sizeMap,
                     mol_setup::MolMap *kindMapFromBox1,
                     mol_setup::SizeMap *sizeMapFromBox1, const uint box);

  static void copyBondInfoIntoMapEntry(const BondAdjacencyList &bondAdjList,
                                       mol_setup::MolMap &kindMap,
                                       std::string fragName);

  // reads BoxTotal PSFs and merges the data, placing the results in kindMap
  // returns 0 if read is successful, -1 on a failure
  int Init(const std::string *psfFilename, const bool *psfDefined,
           pdb_setup::Atoms &pdbAtoms);

  void AssignKinds(const mol_setup::MoleculeVariables &molVars,
                   const FFSetup &ffData);

  // private:
  mol_setup::MolMap kindMap;
  mol_setup::SizeMap sizeMap;
  mol_setup::MoleculeVariables molVars;

private:
  friend class cereal::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &kindMap;
    ar &molVars;
  }
};
#endif
