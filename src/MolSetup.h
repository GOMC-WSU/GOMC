/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOLSETUP_H
#define MOLSETUP_H

#include "BasicTypes.h"
#include "EnsemblePreprocessor.h"

#include <string>
#include <vector>
#include <map>
#include "BondAdjacencyList.h"

namespace config_setup
{
struct RestartSettings;
}
namespace pdb_setup
{
class Atoms;
}
class FFSetup;

namespace mol_setup
{
//!structure to contain an atom's data during initialization
class Atom
{
public:
  Atom(std::string const& l_name, std::string const& l_residue, uint l_resID, std::string const& l_segment, std::string const& l_type,
       const double l_charge, const double l_mass) :
    name(l_name), residue(l_residue), residueID(l_resID), segment(l_segment), type(l_type), charge(l_charge), mass(l_mass) {}
  //private:
  //name (within a molecule) and type (for forcefield params)
  std::string name, type, residue, segment;
  double charge, mass;
  //kind index
  /* ResID is by the PSF Parser to determine multiresidue status by comparing 1st vs all
    of resIDs in a row in the moleculeXAtomIDY 2D vector */
  uint residueID;

  uint kind;
};

class Dihedral
{
public:
  Dihedral(uint atom0, uint atom1, uint atom2, uint atom3)
    : a0(atom0), a1(atom1), a2(atom2), a3(atom3) {}
  //some xplor PSF files have duplicate dihedrals, we need to ignore these
  bool operator == (const Dihedral& other) const;
  bool operator != (const Dihedral& other) const;

  //private:
  //atoms
  uint a0, a1, a2, a3;
  uint kind;
};

class Angle
{
public:
  Angle(uint atom0, uint atom1, uint atom2)
    : a0(atom0), a1(atom1), a2(atom2) {}

  //private:
  uint a0, a1, a2;
  uint kind;
};

class Bond
{
public:
  Bond(uint atom0, uint atom1)
    : a0(atom0), a1(atom1) {}
//   private:
  uint a0, a1;
  uint kind;
};

//!Structure to contain a molecule kind's data during initialization
class MolKind
{
public:
  MolKind() : incomplete(true), isMultiResidue(false), anglesDefined(false), dihedralsDefined(false) {}


  //private:
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;

  uint kindIndex;

  //Used to search PSF file for geometry, meaningless after that
  uint firstAtomID, firstMolID;
  //true while the molecule is still open for modification during PSF read
  bool incomplete;
  bool isMultiResidue;
  bool anglesDefined;
  bool dihedralsDefined;
  std::vector<uint> intraMoleculeResIDs;
};

//List of dihedrals with atom at one end, atom first
std::vector<Dihedral> AtomEndDihs(const MolKind& molKind, uint atom);
//List of dihedrals with atom and partner in middle, atom in a1
std::vector<Dihedral> DihsOnBond(const MolKind& molKind, uint atom, uint partner);
//List of all dihedrals in the molecule kind
std::vector<Dihedral> DihsAll(const MolKind& molKind);
//List of angles with atom at one end, atom first
std::vector<Angle> AtomEndAngles(const MolKind& molKind, uint atom);
//List of angles with atom in middle
std::vector<Angle> AtomMidAngles(const MolKind& molKind, uint atom);
//List of angles with atom at one end, and mid in middle, atom first
std::vector<Angle> AtomMidEndAngles(const MolKind& molKind, uint mid, uint atom);
//List of all angles in the molecule kind
std::vector<Angle> AngsAll(const MolKind& molKind);
//List of bonds with atom at one end, atom first
std::vector<Bond> AtomBonds(const MolKind& molKind, uint atom);
//List of all bonds in the molecule kind
std::vector<Bond> BondsAll(const MolKind& molKind);

//first element (string) is name of molecule type
typedef std::map<std::string, MolKind> MolMap;
typedef std::map<std::size_t, std::vector<std::__cxx11::string> > SizeMap;

/* Used to quickly evaluate if a molecule has been seen by size */
typedef std::map<std::size_t, std::vector<std::__cxx11::string> > SizeMap;
typedef std::map<std::string, uint > KindIndexMap;


//! Reads one or more PSF files into kindMap
/*!
*\param kindMap map to add PSF data to
*\param psfFilename array of strings containing filenames
*\param numFiles number of files to read
*\return -1 if failed, 0 if successful
*/
int ReadCombinePSF(MolMap& kindMap, SizeMap& sizeMap, const std::string* psfFilename,
                   const int numFiles, pdb_setup::Atoms& pdbAtoms);

void PrintMolMapVerbose(const MolMap& kindMap);
void PrintMolMapBrief(const MolMap& kindMap);
}

//wrapper struct for consistent interface
class MolSetup
{
public:
  class Atom;
  int read_atoms(FILE *, unsigned int nAtoms, std::vector<mol_setup::Atom> & allAtoms);
  int createMapAndModifyPDBAtomDataStructure( const BondAdjacencyList & bondAdjList,
                                              const std::vector< std::vector<uint> > & moleculeXAtomIDY, 
                                              std::vector<mol_setup::Atom> & allAtoms,
                                              mol_setup::MolMap & kindMap,
                                              mol_setup::SizeMap & sizeMap,
                                              pdb_setup::Atoms& pdbAtoms,
                                              mol_setup::MolMap * kindMapFromBox1,
                                              mol_setup::SizeMap * sizeMapFromBox1);


  static void copyBondInfoIntoMapEntry(const BondAdjacencyList & bondAdjList, mol_setup::MolMap & kindMap, std::string fragName);


  //reads BoxTotal PSFs and merges the data, placing the results in kindMap
  //returns 0 if read is successful, -1 on a failure
  int Init(const config_setup::RestartSettings& restart,
           const std::string* psfFilename, pdb_setup::Atoms& pdbAtoms);

  void AssignKinds(const pdb_setup::Atoms& pdbAtoms, const FFSetup& ffData);

//private:
  mol_setup::MolMap kindMap;
  mol_setup::SizeMap sizeMap;

};
#endif
