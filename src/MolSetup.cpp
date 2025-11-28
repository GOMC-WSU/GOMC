/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "MolSetup.h"

#include <algorithm> //for swap pre-c++11 compilers
#include <cstdio>
#include <cstring> //strstr
#include <iomanip>
#include <iostream>
#include <sstream> // std::stringstream
#include <utility> //for swap (most modern compilers)

#include "BasicTypes.h"
#include "ConfigSetup.h" //For definition of restart
#include "EnsemblePreprocessor.h"
#include "FFSetup.h" //For geometry kinds
#include "GeomLib.h"
#include "PDBSetup.h" //For mol names->kinds
#include "StrLib.h"

using namespace mol_setup;

bool Dihedral::operator==(const Dihedral &o) const {
  bool same = false;
  if (a0 == o.a0 && a1 == o.a1 && a2 == o.a2 && a3 == o.a3)
    same = true;

  if (a0 == o.a3 && a1 == o.a2 && a2 == o.a1 && a3 == o.a0)
    same = true;

  return same;
}

bool Dihedral::operator!=(const Dihedral &other) const {
  return !(*this == other);
}

bool Improper::operator==(const Improper &o) const {
  bool same = false;
  if (a0 == o.a0 && a1 == o.a1 && a2 == o.a2 && a3 == o.a3)
    same = true;

  // Impropers are order specific, as opposed to Dihedrals
  // if(a0 == o.a3 && a1 == o.a2 && a2 == o.a1 && a3 == o.a0)
  //  same = true;

  return same;
}

bool Improper::operator!=(const Improper &other) const {
  return !(*this == other);
}

namespace {
// Assigns numerical mol kind indices to all molKinds
void AssignMolKinds(MolKind &kind, const mol_setup::MoleculeVariables &molVars,
                    const std::string &name);
void AssignAtomKinds(MolKind &kind, const FFSetup &ffData);
void AssignBondKinds(MolKind &kind, const FFSetup &ffData);
void AssignAngleKinds(MolKind &kind, const FFSetup &ffData);
void AssignDihKinds(MolKind &kind, const FFSetup &ffData);

void BriefBondKinds(MolKind &kind, const FFSetup &ffData);
void BriefAngleKinds(MolKind &kind, const FFSetup &ffData);
void BriefDihKinds(MolKind &kind, const FFSetup &ffData);

// Builds kindMap from PSF file (does not include coordinates) kindMap
// should be empty returns number of atoms in the file, or errors::READ_ERROR if
// the read failed somehow
int ReadPSF(const char *psfFilename, const uint box, MoleculeVariables &molVars,
            MolMap &kindMap, SizeMap &sizeMap, MolMap *kindMapFromBox1 = NULL,
            SizeMap *sizeMapFromBox1 = NULL);
// adds atoms and molecule data in psf to kindMap
// pre: stream is at !NATOMS   post: stream is at end of atom section
int ReadPSFAtoms(FILE *, unsigned int nAtoms,
                 std::vector<mol_setup::Atom> &allAtoms,
                 MoleculeVariables &molVars);
// adds bonds in psf to kindMap
// pre: stream is before !BONDS   post: stream is in bond section just after
// the first appearance of the last molecule
// int ReadPSFBonds(FILE* psf, MolMap& kindMap,
// std::vector<std::pair<uint, std::string> >& firstAtom,
// const uint nbonds);
// adds angles in psf to kindMap
// pre: stream is before !NTHETA   post: stream is in angle section just
// after the first appearance of the last molecule
int ReadPSFAngles(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<uint, std::string>> &firstAtom,
                  const uint nangles);
// adds dihedrals in psf to kindMap
// pre: stream is before !NPHI   post: stream is in dihedral section just
// after the first appearance of the last molecule
int ReadPSFDihedrals(FILE *psf, MolMap &kindMap,
                     std::vector<std::pair<uint, std::string>> &firstAtom,
                     const uint ndihedrals);
// adds impropers in psf to kindMap
// pre: stream is before !NPHI   post: stream is in dihedral section just
// after the first appearance of the last molecule
int ReadPSFImpropers(FILE *psf, MolMap &kindMap,
                     std::vector<std::pair<uint, std::string>> &firstAtom,
                     const uint ndihedrals);
// adds donors in psf to kindMap
// pre: stream is before !NDON   post: stream is in acceptors (NACC) section
// just after the first appearance of the last donor
//
int ReadPSFDonors(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<uint, std::string>> &firstAtom,
                  const uint nDonors);

// adds acceptors in psf to kindMap
// pre: stream is before !NACC   post: stream is in explicit nonbond exclusions
// (NNB) section just after the first appearance of the last donor
//
int ReadPSFAcceptors(FILE *psf, MolMap &kindMap,
                     std::vector<std::pair<uint, std::string>> &firstAtom,
                     const uint nAcceptors);

// adds explicit nonbond exclusions in psf to kindMap
// pre: stream is before !NNB   post: stream is in groups (NGRP) section just
// after the first appearance of the last donor
//
int ReadPSFExplicitNonbondExclusions(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<uint, std::string>> &firstAtom,
    const uint nNonbondExclusions);
// adds groups in psf to kindMap
// pre: stream is before !NGRP   post: stream is in cross-terms (NCRTERM)
// section just after the first appearance of the last donor
//
int ReadPSFGroups(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<uint, std::string>> &firstAtom,
                  const uint nGroups);
// adds cross terms in psf to kindMap
// pre: stream is before !NGRP   post: stream is in cross-terms (NCRTERM)
// section just after the first appearance of the last donor
//
int ReadPSFCrossTerms(FILE *psf, MolMap &kindMap,
                      std::vector<std::pair<uint, std::string>> &firstAtom,
                      const uint nCrossTerms);

} // namespace

// List of dihedrals with atom at one end, atom first
std::vector<Dihedral> mol_setup::AtomEndDihs(const MolKind &molKind,
                                             uint atom) {
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end();
       it < end; ++it) {
    if (it->a0 == atom || it->a3 == atom) {
      result.push_back(*it);
    }
    if (it->a3 == atom) {
      std::swap(result.back().a0, result.back().a3);
      std::swap(result.back().a1, result.back().a2);
    }
  }
  return result;
}

std::vector<Dihedral> mol_setup::DihsOnBond(const MolKind &molKind, uint atom,
                                            uint partner) {
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end();
       it < end; ++it) {
    if (it->a1 == atom && it->a2 == partner) {
      result.push_back(*it);
    } else if (it->a2 == atom && it->a1 == partner) {
      result.push_back(*it);
      std::swap(result.back().a0, result.back().a3);
      std::swap(result.back().a1, result.back().a2);
    }
  }
  return result;
}

std::vector<Dihedral> mol_setup::DihsAll(const MolKind &molKind) {
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end();
       it < end; ++it) {
    result.push_back(*it);
  }
  return result;
}

// List of angles with atom at one end, atom first
std::vector<Angle> mol_setup::AtomEndAngles(const MolKind &molKind, uint atom) {
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end;
       ++it) {
    if (it->a0 == atom || it->a2 == atom) {
      result.push_back(*it);
    }
    if (it->a2 == atom) {
      std::swap(result.back().a0, result.back().a2);
    }
  }
  return result;
}

std::vector<Angle> mol_setup::AtomMidAngles(const MolKind &molKind, uint atom) {
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end;
       ++it) {
    if (it->a1 == atom) {
      result.push_back(*it);
    }
  }
  return result;
}

// List of angles with atom at one end, and mid in middle, atom first
std::vector<Angle> mol_setup::AtomMidEndAngles(const MolKind &molKind, uint mid,
                                               uint atom) {
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end;
       ++it) {
    if ((it->a0 == atom || it->a2 == atom) && (it->a1 == mid)) {
      result.push_back(*it);

      if (it->a2 == atom) {
        std::swap(result.back().a0, result.back().a2);
      }
    }
  }
  return result;
}

std::vector<Angle> mol_setup::AngsAll(const MolKind &molKind) {
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end;
       ++it) {
    result.push_back(*it);
  }
  return result;
}

// List of bonds with atom at one end, atom first
std::vector<Bond> mol_setup::AtomBonds(const MolKind &molKind, uint atom) {
  std::vector<Bond> result;
  typedef std::vector<Bond>::const_iterator Biter;
  for (Biter it = molKind.bonds.begin(), end = molKind.bonds.end(); it < end;
       ++it) {
    if (it->a0 == atom || it->a1 == atom) {
      result.push_back(*it);
    }
    if (it->a1 == atom) {
      std::swap(result.back().a0, result.back().a1);
    }
  }
  return result;
}

std::vector<Bond> mol_setup::BondsAll(const MolKind &molKind) {
  std::vector<Bond> result;
  typedef std::vector<Bond>::const_iterator Biter;
  for (Biter it = molKind.bonds.begin(), end = molKind.bonds.end(); it < end;
       ++it) {
    result.push_back(*it);
  }
  return result;
}

int mol_setup::ReadCombinePSF(MoleculeVariables &molVars, MolMap &kindMap,
                              SizeMap &sizeMap,
                              std::string const *const psfFilename,
                              const bool *psfDefined,
                              pdb_setup::Atoms &pdbAtoms) {
  uint box_0 = 0;
  int errorcode =
      ReadPSF(psfFilename[box_0].c_str(), box_0, molVars, kindMap, sizeMap);
  int nAtoms = errorcode;
  if (errorcode < 0)
    return errorcode;
  MolMap map2;
  SizeMap sizeMap2;
  if ((int)pdbAtoms.count != nAtoms && BOX_TOTAL == 2 && psfDefined[1]) {
    map2.clear();
    uint box_1 = 1;
    errorcode = ReadPSF(psfFilename[box_1].c_str(), box_1, molVars, map2,
                        sizeMap2, &kindMap, &sizeMap);
    nAtoms += errorcode;
    if (errorcode < 0)
      return errorcode;
    kindMap.insert(map2.begin(), map2.end());
  }

  if ((int)pdbAtoms.count != nAtoms) {
    std::cout
        << "Error: This number of atoms in coordinate file(s) (PDB) "
        << pdbAtoms.count
        << " does not match the number of atoms in structure file(s) (PSF) "
        << nAtoms << "!" << std::endl;
    exit(EXIT_FAILURE);
  }

  PrintMolMapVerbose(kindMap);
  // PrintMolMapBrief(kindMap);

  return 0;
}
int MolSetup::Init(const std::string *psfFilename, const bool *psfDefined,
                   pdb_setup::Atoms &pdbAtoms) {
  kindMap.clear();
  sizeMap.clear();
  return ReadCombinePSF(molVars, kindMap, sizeMap, psfFilename, psfDefined,
                        pdbAtoms);
}

void MolSetup::AssignKinds(const mol_setup::MoleculeVariables &molVars,
                           const FFSetup &ffData) {
  typedef MolMap::iterator MapIt;
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    /* Indices determined using molVars name vector,
    so we pass molVars.moleculeKindNames to outputVars
    in initialization of CPUSide */
    AssignMolKinds(it->second, molVars, it->first);
    AssignAtomKinds(it->second, ffData);
    AssignBondKinds(it->second, ffData);
    AssignAngleKinds(it->second, ffData);
    AssignDihKinds(it->second, ffData);
  }

  // Print bonded Information
  printf("Bonds parameter:\n");
  printf("%s %33s %15s \n", "Atom Types", "Kb(K)", "b0(A)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefBondKinds(it->second, ffData);
  }

  printf("Angles parameter:\n");
  printf("%s %33s %22s \n", "Atom Types", "Ktheta(K)", "theta0(degree)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefAngleKinds(it->second, ffData);
  }

  printf("Dihedrals parameter:\n");
  printf("%s %33s %4s %16s \n", "Atom Types", "Kchi(K)", "n", "delta(degree)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefDihKinds(it->second, ffData);
  }
  std::cout << std::endl;
}

typedef std::vector<uint>::const_iterator candidateIterator;

void createKindMap(mol_setup::MoleculeVariables &molVars,
                   const BondAdjacencyList &bondAdjList,
                   const std::vector<std::vector<uint>> &moleculeXAtomIDY,
                   std::vector<mol_setup::Atom> &allAtoms,
                   mol_setup::MolMap &kindMap, mol_setup::SizeMap &sizeMap,
                   mol_setup::MolMap *kindMapFromBox1,
                   mol_setup::SizeMap *sizeMapFromBox1, const uint box) {
  uint startIdxAtomBoxOffset;
  AlphaNum uniqueSuffixGenerator;
  if (box == 0) {
    startIdxAtomBoxOffset = 0;
  } else {
    if (molVars.numberMolsInBox0 != 0) {
      startIdxAtomBoxOffset = molVars.lastAtomIndexInBox0 + 1;
    } else {
      startIdxAtomBoxOffset = 0;
    }
  }

  for (std::vector<std::vector<uint>>::const_iterator it =
           moleculeXAtomIDY.cbegin();
       it != moleculeXAtomIDY.cend(); it++) {
    std::string fragName;
    bool foundEntryInOldMap = false;

    if (sizeMapFromBox1 != NULL && kindMapFromBox1 != NULL) {
      /* A size -> moleculeKind map for quick evaluation of new molecules based
        on molMap entries of a given size exisitng or not */
      /* Search by size for existing molecules from Box 1 if it exists*/
      SizeMap::iterator sizeIt = sizeMapFromBox1->find(it->size());

      // Found a match in the old molMap.
      if (sizeIt != sizeMapFromBox1->end()) {
        /* Iterate through all the size consistent map entries */
        for (std::vector<std::string>::const_iterator sizeConsistentEntries =
                 sizeIt->second.cbegin();
             sizeConsistentEntries != sizeIt->second.cend();
             sizeConsistentEntries++) {
          /* Iterate atom by atom of a given size consistent map entries with
           * the candidate molecule*/
          typedef std::vector<mol_setup::Atom>::const_iterator atomIterator;
          std::pair<atomIterator, candidateIterator> itPair(
              (*kindMapFromBox1)[*sizeConsistentEntries].atoms.cbegin(),
              it->cbegin());
          for (; itPair.second != it->cend(); ++itPair.first, ++itPair.second) {
            /* Atom equality operator
              Checks atom name, type, residue, charge, and mass */
            if (*(itPair.first) == allAtoms[*(itPair.second)]) {
              continue;
            } else {
              break;
            }
          }

          /* Found no matching molecules in Box 2 map */
          if (itPair.second == it->cend()) {
            /* Get the map key */
            fragName = *sizeConsistentEntries;
            /* Boilerplate PDB Data modifications for matches */
            molVars.startIdxMolecules.push_back(startIdxAtomBoxOffset +
                                                it->front());
            molVars.moleculeKinds.push_back(
                (*kindMapFromBox1)[fragName].kindIndex);
            molVars.moleculeNames.push_back(
                (*kindMapFromBox1)[fragName].isMultiResidue
                    ? fragName
                    : ((*kindMapFromBox1)[fragName].moleculeName));
            molVars.moleculeSegmentNames.push_back(
                allAtoms[it->front()].segment);

            /* Boilerplate PDB Data modifications for matches */

            /* Search current KindMap for this entry.
              We won't have a value for fragName unless we matched
              in the old molMap, since it may be a PROTX name.
              That's why this all has to be in this conditional.
            */
            MolMap::const_iterator kindIt = kindMap.find(fragName);

            /* If we don't have it in our new Map, then we need to add the first
              entry in our new Map to quench the post processessing
              requirements. */
            if (kindIt == kindMap.cend()) {
              if ((*kindMapFromBox1)[fragName].isMultiResidue) {
                kindMap[fragName] = MolKind();
                kindMap[fragName].isMultiResidue = true;
                uint intraResID = 0;
                uint compareResID = allAtoms[it->front()].residueID;
                for (std::vector<uint>::const_iterator connectedComponentIt =
                         it->cbegin();
                     connectedComponentIt != it->cend();
                     connectedComponentIt++) {
                  kindMap[fragName].atoms.push_back(
                      allAtoms[*connectedComponentIt]);
                  if (compareResID ==
                      allAtoms[*connectedComponentIt].residueID) {
                    kindMap[fragName].intraMoleculeResIDs.push_back(intraResID);
                  } else {
                    compareResID = allAtoms[*connectedComponentIt].residueID;
                    ++intraResID;
                    kindMap[fragName].intraMoleculeResIDs.push_back(intraResID);
                  }
                }
              } else {
                kindMap[fragName] = MolKind();
                kindMap[fragName].isMultiResidue = false;
                for (std::vector<uint>::const_iterator connectedComponentIt =
                         it->cbegin();
                     connectedComponentIt != it->cend();
                     connectedComponentIt++) {
                  kindMap[fragName].atoms.push_back(
                      allAtoms[*connectedComponentIt]);
                }
              }
              kindMap[fragName].firstAtomID = it->front() + 1;
              kindMap[fragName].firstMolID = allAtoms[it->front()].residueID;
              kindMap[fragName].kindIndex =
                  (*kindMapFromBox1)[fragName].kindIndex;
              kindMap[fragName].moleculeName =
                  (*kindMapFromBox1)[fragName].isMultiResidue
                      ? fragName
                      : (*kindMapFromBox1)[fragName].moleculeName;
              MolSetup::copyBondInfoIntoMapEntry(bondAdjList, kindMap,
                                                 fragName);

              /* Finally, search new sizeMap for existing entry of same size as
                 old molecule. This handles chance that we have two equal sized,
                 but different, molecules in
                 the two boxes.  Since the map entry will already exist, it
                 wouldn't be a newSize */

              SizeMap::iterator sizeIt = sizeMap.find(it->size());
              /* New Size */
              if (sizeIt == sizeMap.end()) {
                sizeMap[it->size()] = std::vector<std::string>{fragName};
              } else {
                sizeMap[it->size()].push_back(fragName);
              }
            }
            foundEntryInOldMap = true;
            /* We found our entry, so there is no need to continue iterating
            through size consistent entries */
            break;
          }
        }
      }
    }

    /* Follow procedure from earlier versions of single box ensemble map
       entries. Very similar to above, except we may need to generate fragName
       in the case of a protein. */
    if (!foundEntryInOldMap) {
      /* Search by size for existing molecules  for quick evaluation of new
      molecules based on molMap entries of a given size exisitng or not */
      SizeMap::iterator sizeIt = sizeMap.find(it->size());
      std::string fragName;
      bool multiResidue = false;
      bool newSize = false;
      bool newMapEntry = true;

      typedef std::vector<uint>::const_iterator candidateIterator;
      /* Found no matching molecules by size */
      if (sizeIt == sizeMap.end()) {
        newSize = true;
        /* Found molecules with the same size, now evaluate for atom equality */
      } else {
        /* Iterate through all the size consistent map entries */
        for (std::vector<std::string>::const_iterator sizeConsistentEntries =
                 sizeIt->second.cbegin();
             sizeConsistentEntries != sizeIt->second.cend();
             sizeConsistentEntries++) {
          /* Iterate atom by atom of a given size consistent map entries with
           * the candidate molecule*/
          typedef std::vector<mol_setup::Atom>::const_iterator atomIterator;
          std::pair<atomIterator, candidateIterator> itPair(
              kindMap[*sizeConsistentEntries].atoms.cbegin(), it->cbegin());
          for (; itPair.second != it->cend(); ++itPair.first, ++itPair.second) {
            if (*(itPair.first) == allAtoms[*(itPair.second)]) {
              continue;
            } else {
              break;
            }
          }
          // Found a match
          if (itPair.second == it->cend()) {
            fragName = *sizeConsistentEntries;
            // Modify PDBData
            molVars.startIdxMolecules.push_back(startIdxAtomBoxOffset +
                                                it->front());
            molVars.moleculeKinds.push_back(kindMap[fragName].kindIndex);
            molVars.moleculeNames.push_back(
                kindMap[fragName].isMultiResidue
                    ? fragName
                    : (kindMap[fragName].moleculeName));
            molVars.moleculeSegmentNames.push_back(
                allAtoms[it->front()].segment);

            newMapEntry = false;
            break;
          }
        }
      }

      if (newMapEntry) {
        /* Determine if this connected component is a standard or multiResidue
         * molecule */
        for (candidateIterator connectedComponentIt = it->cbegin();
             connectedComponentIt != it->cend(); connectedComponentIt++) {
          if (allAtoms[*connectedComponentIt].residueID ==
              allAtoms[it->front()].residueID) {
            continue;
          } else {
            multiResidue = true;
            break;
          }
        }
        if (multiResidue) {
          fragName = "PROT" + uniqueSuffixGenerator.uint2String(
                                  molVars.stringSuffixMultiResidue);
          molVars.stringSuffixMultiResidue++;
          printf("\n%-40s \n",
                 "Warning: A molecule containing > 1 residue is detected.");
          printf("The simulation will name it %s.\n", fragName.c_str());
          printf("See the chart at the end of the output log describing this "
                 "entry.\n");

          kindMap[fragName] = MolKind();
          kindMap[fragName].isMultiResidue = true;
          kindMap[fragName].moleculeName = fragName;
          uint intraResID = 0;
          uint compareResID = allAtoms[it->front()].residueID;
          for (std::vector<uint>::const_iterator connectedComponentIt =
                   it->cbegin();
               connectedComponentIt != it->cend(); connectedComponentIt++) {
            kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
            if (compareResID == allAtoms[*connectedComponentIt].residueID) {
              kindMap[fragName].intraMoleculeResIDs.push_back(intraResID);
            } else {
              compareResID = allAtoms[*connectedComponentIt].residueID;
              ++intraResID;
              kindMap[fragName].intraMoleculeResIDs.push_back(intraResID);
            }
          }
        } else {
          // Generate Unique MapKey instead of using residue
          // fragName = allAtoms[it->front()].residue;
          fragName = "MOL" + uniqueSuffixGenerator.uint2String(
                                 molVars.stringSuffixNonMultiResidue);
          molVars.stringSuffixNonMultiResidue++;
          kindMap[fragName] = MolKind();
          kindMap[fragName].isMultiResidue = false;
          kindMap[fragName].moleculeName = allAtoms[it->front()].residue;
          for (std::vector<uint>::const_iterator connectedComponentIt =
                   it->cbegin();
               connectedComponentIt != it->cend(); connectedComponentIt++) {
            kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
          }
        }
        kindMap[fragName].firstAtomID = it->front() + 1;
        kindMap[fragName].firstMolID = allAtoms[it->front()].residueID;
        kindMap[fragName].kindIndex = molVars.molKindIndex;
        molVars.startIdxMolecules.push_back(startIdxAtomBoxOffset +
                                            kindMap[fragName].firstAtomID - 1);
        molVars.moleculeKinds.push_back(kindMap[fragName].kindIndex);
        molVars.moleculeKindNames.push_back(
            kindMap[fragName].isMultiResidue
                ? fragName
                : (kindMap[fragName].moleculeName));
        molVars.uniqueMapKeys.push_back(fragName);
        molVars.moleculeNames.push_back(kindMap[fragName].isMultiResidue
                                            ? fragName
                                            : (kindMap[fragName].moleculeName));
        molVars.moleculeSegmentNames.push_back(allAtoms[it->front()].segment);
        MolSetup::copyBondInfoIntoMapEntry(bondAdjList, kindMap, fragName);
        molVars.molKindIndex++;
        if (newSize) {
          sizeMap[it->size()] = std::vector<std::string>{fragName};
        } else {
          sizeMap[it->size()].push_back(fragName);
        }
      }
    }
    molVars.moleculeIteration++;
  }
  if (box == 0) {
    molVars.numberMolsInBox0 = moleculeXAtomIDY.size();
    if (molVars.numberMolsInBox0 != 0)
      molVars.lastAtomIndexInBox0 = (moleculeXAtomIDY.back()).back();
  }
}

typedef std::map<std::string, mol_setup::MolKind> MolMap;
void MolSetup::copyBondInfoIntoMapEntry(const BondAdjacencyList &bondAdjList,
                                        mol_setup::MolMap &kindMap,
                                        std::string fragName) {
  int molBegin = static_cast<int>(kindMap[fragName].firstAtomID) - 1;
  // index AFTER last atom in molecule
  int molEnd = molBegin + static_cast<int>(kindMap[fragName].atoms.size());
  // assign the bond
  for (int i = molBegin; i < molEnd; i++) {
    adjNode *ptr = bondAdjList.head[i];
    while (ptr != nullptr) {
      if (i < ptr->val) {
        kindMap[fragName].bonds.push_back(
            Bond(i - molBegin, ptr->val - molBegin));
      }
      ptr = ptr->next;
    }
  }
  /* before returning, reverse the bonds vector, which should get dead on balls
   * accurate merged psfs across cycles */
  std::reverse(kindMap[fragName].bonds.begin(), kindMap[fragName].bonds.end());
}

namespace {

void AssignMolKinds(MolKind &kind, const mol_setup::MoleculeVariables &molVars,
                    const std::string &name) {
  uint index = std::find(molVars.moleculeKindNames.begin(),
                         molVars.moleculeKindNames.end(), name) -
               molVars.moleculeKindNames.begin();
  kind.kindIndex = index;
}

void AssignAtomKinds(MolKind &kind, const FFSetup &ffData) {
  for (uint i = 0; i < kind.atoms.size(); ++i) {
    int thisKind = ffData.mie.Find(&kind.atoms[i].type, ffData.mie.name);
    if (thisKind < 0) {
      fprintf(stderr,
              "ERROR: Atom Type %s not specified in nonbonded section of "
              "parameter file.\n",
              kind.atoms[i].type.c_str());
      exit(EXIT_FAILURE);
    }
    kind.atoms[i].kind = thisKind;
  }
}

void AssignBondKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 2;
  std::string elementNames[ATOMS_PER];

  int search = 0;
  for (uint i = 0; i < kind.bonds.size(); ++i) {
    elementNames[0] = kind.atoms[kind.bonds[i].a0].type;
    elementNames[1] = kind.atoms[kind.bonds[i].a1].type;
    search = ffData.bond.Find(elementNames, ffData.bond.name);
    if (search >= 0) {
      kind.bonds[i].kind = search;
    } else {
      std::string missing;
      for (uint m = 0; m < ATOMS_PER; ++m)
        missing.append(elementNames[m]).append(" ");
      fprintf(stderr, "Error: Bond %s not found in parameter files.\n",
              missing.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

void BriefBondKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 2;
  std::string elementNames[ATOMS_PER];
  std::vector<std::string> printed;

  if (kind.bonds.size() == 0)
    return;

  for (uint i = 0; i < kind.bonds.size(); ++i) {
    uint search = kind.bonds[i].kind;
    std::string bondName;

    elementNames[0] = kind.atoms[kind.bonds[i].a0].type;
    elementNames[1] = kind.atoms[kind.bonds[i].a1].type;

    for (uint m = 0; m < ATOMS_PER; ++m) {
      bondName.append(elementNames[m]).append("\t");
    }

    if (find(printed.begin(), printed.end(), bondName) == printed.end()) {
      printf("%s", bondName.c_str());
      if (ffData.bond.GetKb(search) > 99999999)
        printf("%28s %16.4f \n", "FIX", ffData.bond.Getb0(search));
      else
        printf("%28.4f %16.4f \n", ffData.bond.GetKb(search),
               ffData.bond.Getb0(search));

      printed.push_back(bondName);
    }
  }
  std::cout << std::endl;
}

void AssignAngleKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 3;
  std::string elementNames[ATOMS_PER];

  int search = 0;
  for (uint i = 0; i < kind.angles.size(); ++i) {
    elementNames[0] = kind.atoms[kind.angles[i].a0].type;
    elementNames[1] = kind.atoms[kind.angles[i].a1].type;
    elementNames[2] = kind.atoms[kind.angles[i].a2].type;
    search = ffData.angle.Find(elementNames, ffData.angle.name);
    if (search >= 0) {
      kind.angles[i].kind = search;
    } else {
      std::string missing;
      for (uint m = 0; m < ATOMS_PER; ++m)
        missing.append(elementNames[m]).append(" ");
      fprintf(stderr, "Error: Angle %s not found in parameter files.\n",
              missing.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

void BriefAngleKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 3;
  std::string elementNames[ATOMS_PER];
  std::vector<std::string> printed;
  double coef = 180.0 * M_1_PI;

  if (kind.angles.size() == 0)
    return;

  for (uint i = 0; i < kind.angles.size(); ++i) {
    std::string angleName;
    uint search = kind.angles[i].kind;
    elementNames[0] = kind.atoms[kind.angles[i].a0].type;
    elementNames[1] = kind.atoms[kind.angles[i].a1].type;
    elementNames[2] = kind.atoms[kind.angles[i].a2].type;

    for (uint m = 0; m < ATOMS_PER; ++m) {
      angleName.append(elementNames[m]).append("\t");
    }

    if (find(printed.begin(), printed.end(), angleName) == printed.end()) {
      printf("%s", angleName.c_str());
      if (ffData.angle.GetKtheta(search) > 99999999)
        printf("%20s %16.4f \n", "FIX", ffData.angle.Gettheta0(search) * coef);
      else
        printf("%20.4f %16.4f \n", ffData.angle.GetKtheta(search),
               ffData.angle.Gettheta0(search) * coef);

      printed.push_back(angleName);
    }
  }
  std::cout << std::endl;
}

void AssignDihKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 4;
  std::string elementNames[ATOMS_PER];

  int search = 0;
  for (uint i = 0; i < kind.dihedrals.size(); ++i) {
    elementNames[0] = kind.atoms[kind.dihedrals[i].a0].type;
    elementNames[1] = kind.atoms[kind.dihedrals[i].a1].type;
    elementNames[2] = kind.atoms[kind.dihedrals[i].a2].type;
    elementNames[3] = kind.atoms[kind.dihedrals[i].a3].type;
    search = ffData.dih.Find(elementNames, ffData.dih.name);
    if (search < 0) {
      std::string missing;
      for (uint m = 0; m < ATOMS_PER; ++m)
        missing.append(elementNames[m]).append(" ");
      fprintf(stderr, "Error: Dihedral %s not found in parameter files.\n",
              missing.c_str());
      exit(EXIT_FAILURE);
    }
    kind.dihedrals[i].kind = search;
  }
}

void BriefDihKinds(MolKind &kind, const FFSetup &ffData) {
  const uint ATOMS_PER = 4;
  std::string elementNames[ATOMS_PER];
  double coef = 180.0 * M_1_PI;
  std::vector<std::string> printed;

  if (kind.dihedrals.size() == 0)
    return;

  for (uint i = 0; i < kind.dihedrals.size(); ++i) {
    std::string dName = ffData.dih.name[kind.dihedrals[i].kind];
    std::string dihedralName;
    uint dihsize = ffData.dih.GetSizeDih(dName);

    elementNames[0] = kind.atoms[kind.dihedrals[i].a0].type;
    elementNames[1] = kind.atoms[kind.dihedrals[i].a1].type;
    elementNames[2] = kind.atoms[kind.dihedrals[i].a2].type;
    elementNames[3] = kind.atoms[kind.dihedrals[i].a3].type;

    for (uint m = 0; m < ATOMS_PER; ++m) {
      dihedralName.append(elementNames[m]).append("\t");
    }

    if (find(printed.begin(), printed.end(), dihedralName) == printed.end()) {
      for (uint j = 0; j < dihsize; j++) {
        printf("%s", dihedralName.c_str());
        printf("%12.4f %4d %11.4f \n", ffData.dih.GetKchi(dName, j),
               ffData.dih.Getn(dName, j), ffData.dih.Getdelta(dName, j) * coef);
      }
      printed.push_back(dihedralName);
    }
  }
}

} // namespace

void mol_setup::PrintMolMapVerbose(const MolMap &kindMap) {
  std::cout << "\nMolecules in PSF:\n";
  MolMap::const_iterator it = kindMap.begin();
  while (it != kindMap.end()) {
    std::cout << "Molecule Kind: " << it->second.moleculeName << std::endl;
    std::cout << "Idx\tname\ttype\tcharge\tmass\n";
    for (uint i = 0; i < it->second.atoms.size(); i++) {
      std::cout << i << "\t" << it->second.atoms[i].name << '\t'
                << it->second.atoms[i].type << '\t' << std::setprecision(4)
                << it->second.atoms[i].charge << '\t' << std::setprecision(4)
                << it->second.atoms[i].mass << std::endl;
    }
    std::cout << "\nBonds:";
    for (uint i = 0; i < it->second.bonds.size(); i++) {
      if (i % 10 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.bonds[i].a0 << ' '
                << it->second.bonds[i].a1 << ']' << ' ';
    }
    std::cout << std::endl << "\nAngles:";
    for (uint i = 0; i < it->second.angles.size(); i++) {
      if (i % 7 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.angles[i].a0 << ' '
                << it->second.angles[i].a1 << ' ' << it->second.angles[i].a2
                << ']' << ' ';
    }
    std::cout << std::endl << "\nDihedrals:";
    for (uint i = 0; i < it->second.dihedrals.size(); i++) {
      if (i % 5 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.dihedrals[i].a0 << ' '
                << it->second.dihedrals[i].a1 << ' '
                << it->second.dihedrals[i].a2 << ' '
                << it->second.dihedrals[i].a3 << ']' << ' ';
    }
    std::cout << std::endl << "\nImpropers:";
    for (uint i = 0; i < it->second.impropers.size(); i++) {
      if (i % 5 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.impropers[i].a0 << ' '
                << it->second.impropers[i].a1 << ' '
                << it->second.impropers[i].a2 << ' '
                << it->second.impropers[i].a3 << ']' << ' ';
    }
    ++it;
    std::cout << std::endl << std::endl;
  }
}

void mol_setup::PrintMolMapBrief(const MolMap &kindMap) {
  std::cout << "Molecules in PSF:\n";
  std::cout << "Name\t#Atom\t#Bond\t#Ang\t#Dih\t\n";
  MolMap::const_iterator it = kindMap.begin();
  while (it != kindMap.end()) {
    std::cout << it->first << '\t' << it->second.atoms.size() << '\t'
              << it->second.bonds.size() << '\t' << it->second.angles.size()
              << '\t' << it->second.dihedrals.size() << '\t'
              << it->second.impropers.size() << '\t' << std::endl;
    ++it;
  }
}

namespace {
// Initializes system from PSF file (does not include coordinates)
// returns number of atoms in the file, or errors::READ_ERROR if the read failed
// somehow
int ReadPSF(const char *psfFilename, const uint box, MoleculeVariables &molVars,
            MolMap &kindMap, SizeMap &sizeMap, MolMap *kindMapFromBox1,
            SizeMap *sizeMapFromBox1) {
  FILE *psf = fopen(psfFilename, "r");
  char *check; // return value of fgets
  int count;   // for number of bonds/angles/dihs
  if (psf == NULL) {
    fprintf(
        stderr,
        "ERROR: Failed to open PSF file %s for molecule data.\nExiting...\n",
        psfFilename);
    return errors::READ_ERROR;
  }
  char input[512];
  unsigned int nAtoms;
  // find atom header+count
  do {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read atoms from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  } while (strstr(input, "!NATOM") == NULL);
  sscanf(input, " %u", &nAtoms);

  /* GJS - This is a flat vector of atom objects.  Instead of building kindMap
  entries as we parse the PSF file, we will first read the entire atom section,
  building an N length array of atoms, N is number of atoms.  Since, we process
  the PSF file as a flat file handle, we cannot process bonds before atoms
  without physically generating a new PSFFile reversing the order of ATOMS <->
  BONDS.
  Hence, the necessity to build this vector before knowing how the atoms are
  connected. */
  std::vector<mol_setup::Atom> allAtoms;
  ReadPSFAtoms(psf, nAtoms, allAtoms, molVars);
  // build list of start particles for each type, so we can find it and skip
  // everything else
  // make sure molecule has bonds, appears before !NBOND
  // find bond header+count
  while (strstr(input, "!NBOND") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read bonds from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  count = atoi(input);
  /* moleculeXAtomIDY - A 2D vector, of jagged row length,
    of atom indices into the allAtoms vector, with atoms indexed from 0.
    All indices in a row will make up a connected component.  Note, this is
    a redundant vector, and there may be repeated molecules.  These will be
    consolidated and counted in the (1)createMap portion of
    createMapAndModifyPDBAtomDataStructure.

            Y - atom 1 atom 2 .. atom N
  X - molecule 1  0     3
  .
  .
  .
  X - molecule M
  */
  std::vector<std::vector<uint>> moleculeXAtomIDY;
  /* A standard adjacency list with N nodes, where N is number of atoms.
   This is an undirected graph, where edges between nodes
   represent bonds between atoms.  It is generated by DFS, checking if a node
   has been visited before.
  */
  BondAdjacencyList bondAdjList(psf, nAtoms, count, moleculeXAtomIDY);
  /* Molecular equality is determined by a series of evaluations, with early
    termination - 1) Length of candidate versus all entries in the map 2) Atom
    by atom check for equality between candidate and all entries.

    If the candidate matches an exisiting map entry, we push the first resID
    onto startIDxRes, along with other necessary information.

    If the candidate is determined to be novel, we first determine its status as
    standard or multiresidue.

    If multiresidue a dummy string is generated and atoms added to the entry
    along with intraMolecularResIDs normalized to 0 .. n-1, n is number of
    residues in the multiresidue molecule.

    Otherwise, standard procedure for creating a map entry is followed.

    The bond information contained in the Adjacency list is assigned to map
    entries.

  */
  createKindMap(molVars, bondAdjList, moleculeXAtomIDY, allAtoms, kindMap,
                sizeMap, kindMapFromBox1, sizeMapFromBox1, box);

  std::vector<std::pair<unsigned int, std::string>> firstAtomLookup;
  for (MolMap::iterator it = kindMap.begin(); it != kindMap.end(); ++it) {
    firstAtomLookup.push_back(
        std::make_pair(it->second.firstAtomID, it->first));
  }
  std::sort(firstAtomLookup.begin(), firstAtomLookup.end());
  // find bond header+count
  // make sure molecule has bonds, appears before !NBOND
  // find angle header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NTHETA") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read angles from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  // make sure molecule has angles, count appears before !NTHETA
  count = atoi(input);
  if (ReadPSFAngles(psf, kindMap, firstAtomLookup, count) ==
      errors::READ_ERROR) {
    fclose(psf);
    return errors::READ_ERROR;
  }
  // find dihedrals header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NPHI") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read dihedrals from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  // make sure molecule has dihs, count appears before !NPHI
  count = atoi(input);
  if (ReadPSFDihedrals(psf, kindMap, firstAtomLookup, count) ==
      errors::READ_ERROR) {
    fclose(psf);
    return errors::READ_ERROR;
  }

  // find impropers header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NIMPHI") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read impropers from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  // make sure molecule has imps, count appears before !NIMPHI
  count = atoi(input);
  if (ReadPSFImpropers(psf, kindMap, firstAtomLookup, count) ==
      errors::READ_ERROR) {
    fclose(psf);
    return errors::READ_ERROR;
  }

  // find donors header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NDON") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read donors from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  // make sure molecule has donors, count appears before !NDON
  count = atoi(input);
  if (ReadPSFDonors(psf, kindMap, firstAtomLookup, count) ==
      errors::READ_ERROR) {
    fclose(psf);
    return errors::READ_ERROR;
  }

  // find acceptors header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NACC") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read acceptors from PSF file %s\n",
              psfFilename);
      fclose(psf);
      return errors::READ_ERROR;
    }
  }
  // make sure molecule has acceptors, count appears before !NACC
  count = atoi(input);
  if (ReadPSFAcceptors(psf, kindMap, firstAtomLookup, count) ==
      errors::READ_ERROR) {
    fclose(psf);
    return errors::READ_ERROR;
  }
  /*

   //find explicit nonbond exclusions  header+count
   fseek(psf, 0, SEEK_SET);
   while (strstr(input, "!NNB") == NULL) {
     check = fgets(input, 511, psf);
     if (check == NULL) {
       fprintf(stderr, "ERROR: Unable to read explicit nonbond exclusions from
   PSF file %s\n", psfFilename); fclose(psf); return errors::READ_ERROR;
     }
   }
   //make sure molecule has explicit nonbond exclusions, count appears before
   !NNB count = atoi(input); if (ReadPSFExplicitNonbondExclusions(psf, kindMap,
   firstAtomLookup, count) == errors::READ_ERROR) { fclose(psf); return
   errors::READ_ERROR;
   }

   //find groups header+count
   fseek(psf, 0, SEEK_SET);
   while (strstr(input, "!NGRP") == NULL) {
     check = fgets(input, 511, psf);
     if (check == NULL) {
       fprintf(stderr, "ERROR: Unable to read groups from PSF file %s\n",
               psfFilename);
       fclose(psf);
       return errors::READ_ERROR;
     }
   }
   //make sure molecule has groups, count appears before !NGRP
   count = atoi(input);
   if (ReadPSFGroups(psf, kindMap, firstAtomLookup, count) ==
   errors::READ_ERROR) { fclose(psf); return errors::READ_ERROR;
   }

   //find cross terms header+count
   fseek(psf, 0, SEEK_SET);
   while (strstr(input, "!NCRTERM") == NULL) {
     check = fgets(input, 511, psf);
     if (check == NULL) {
       fprintf(stderr, "ERROR: Unable to read cross terms from PSF file %s\n",
               psfFilename);
       fclose(psf);
       return errors::READ_ERROR;
     }
   }
   //make sure molecule has cross terms, count appears before !NCRTERM
   count = atoi(input);
   if (ReadPSFCrossTerms(psf, kindMap, firstAtomLookup, count) ==
   errors::READ_ERROR) { fclose(psf); return errors::READ_ERROR;
   }
 */
  fclose(psf);

  return nAtoms;
}

// adds atoms and molecule data in psf to kindMap
// pre: stream is at !NATOMS   post: stream is at end of atom section
int ReadPSFAtoms(FILE *psf, unsigned int nAtoms,
                 std::vector<mol_setup::Atom> &allAtoms,
                 MoleculeVariables &molVars) {
  char input[512];
  unsigned int atomID = 0;
  unsigned int molID;
  char segment[11], moleculeName[11], atomName[11], atomType[11];
  double charge, mass;

  while (atomID < nAtoms) {
    char *check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Could not find all atoms in PSF file ");
      return errors::READ_ERROR;
    }
    // skip comment/blank lines
    if (input[0] == '!' || str::AllWS(input))
      continue;
    // parse line
    sscanf(input, " %u %s %u %s %s %s %lf %lf ", &atomID, segment, &molID,
           moleculeName, atomName, atomType, &charge, &mass);
    allAtoms.push_back(mol_setup::Atom(atomName, moleculeName, molID, segment,
                                       atomType, charge, mass));
  }
  return 0;
}

// adds angles in psf to kindMap
// pre: stream is before !NTHETA   post: stream is in angle section just after
// the first appearance of the last molecule
int ReadPSFAngles(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<unsigned int, std::string>> &firstAtom,
                  const uint nangles) {
  unsigned int atom0, atom1, atom2;
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nangles; n++) {
    dummy = fscanf(psf, "%u %u %u", &atom0, &atom1, &atom2);
    if (dummy != 3) {
      fprintf(stderr, "ERROR: Incorrect Number of angles in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all angles in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this angles
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign the angle
      if (atom0 >= molBegin && atom0 < molEnd) {
        currentMol.angles.push_back(
            Angle(atom0 - molBegin, atom1 - molBegin, atom2 - molBegin));
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  // Check if we defined all angles
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind &currentMol = kindMap[firstAtom[i].second];
    if (currentMol.atoms.size() > 2 && !defined[i]) {
      std::cout << "Warning: Angle is missing for " << firstAtom[i].second
                << " !\n";
    }
  }
  return 0;
}

// bool ContainsDihedral(const std::vector<uint>& vec, const int dih[])
// {
// for (uint i = 0; i < vec.size(); i += 4) {
// bool match = true;
// for (uint j = 0; j < 4; ++j)
// if ((int) vec[i + j] != dih[j])
// match = false;
// if (match)
// return true;
// }
// return false;
// }

// adds dihedrals in psf to kindMap
// pre: stream is before !NPHI   post: stream is in dihedral section just after
// the first appearance of the last molecule
//
int ReadPSFDihedrals(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<unsigned int, std::string>> &firstAtom,
    const uint ndihedrals) {
  Dihedral dih(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < ndihedrals; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &dih.a0, &dih.a1, &dih.a2, &dih.a3);
    if (dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of dihedrals in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all dihedrals in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this dihedral
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign dihedral
      if (dih.a0 >= molBegin && dih.a0 < molEnd) {
        dih.a0 -= molBegin;
        dih.a1 -= molBegin;
        dih.a2 -= molBegin;
        dih.a3 -= molBegin;
        // some xplor PSF files have duplicate dihedrals, we need to ignore
        // these
        if (std::find(currentMol.dihedrals.begin(), currentMol.dihedrals.end(),
                      dih) == currentMol.dihedrals.end()) {
          currentMol.dihedrals.push_back(dih);
        }
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  // Check if we defined all dihedrals
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind &currentMol = kindMap[firstAtom[i].second];
    if (currentMol.atoms.size() > 3 && !defined[i]) {
      std::cout << "Warning: Dihedral is missing for " << firstAtom[i].second
                << " !\n";
    }
  }
  return 0;
}

// adds impropers in psf to kindMap
// pre: stream is before !NIMPHI   post: stream is in donors (NDON) section just
// after the first appearance of the last improper
//
int ReadPSFImpropers(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<unsigned int, std::string>> &firstAtom,
    const uint nimpropers) {
  Improper imp(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nimpropers; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &imp.a0, &imp.a1, &imp.a2, &imp.a3);
    if (dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of impropers in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all impropers in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this impropers
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign impropers
      if (imp.a0 >= molBegin && imp.a0 < molEnd) {
        imp.a0 -= molBegin;
        imp.a1 -= molBegin;
        imp.a2 -= molBegin;
        imp.a3 -= molBegin;
        // some xplor PSF files have duplicate impropers, we need to ignore
        // these
        if (std::find(currentMol.impropers.begin(), currentMol.impropers.end(),
                      imp) == currentMol.impropers.end()) {
          currentMol.impropers.push_back(imp);
        }
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }

  return 0;
}

// adds donors in psf to kindMap
// pre: stream is before !NDON   post: stream is in acceptors (NACC) section
// just after the first appearance of the last donor
//
int ReadPSFDonors(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<unsigned int, std::string>> &firstAtom,
                  const uint nDonors) {
  unsigned int atom0, atom1;
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nDonors; n++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if (dummy != 2) {
      fprintf(stderr, "ERROR: Incorrect Number of donors in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all donors in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this bond
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in molecule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign the bond
      if (atom0 >= molBegin && atom0 < molEnd) {
        currentMol.donors.push_back(Bond(atom0 - molBegin, atom1 - molBegin));
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }

  return 0;
}

// adds acceptors in psf to kindMap
// pre: stream is before !NACC   post: stream is in explicit nonbond exclusions
// (NNB) section just after the first appearance of the last donor
//
int ReadPSFAcceptors(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<unsigned int, std::string>> &firstAtom,
    const uint nAcceptors) {
  unsigned int atom0, atom1;
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nAcceptors; n++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if (dummy != 2) {
      fprintf(stderr, "ERROR: Incorrect Number of acceptors in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all acceptors in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this bond
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in molecule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign the bond
      if (atom0 >= molBegin && atom0 < molEnd) {
        currentMol.acceptors.push_back(
            Bond(atom0 - molBegin, atom1 - molBegin));
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }

  return 0;
}

/* Explanation of NNB Format

Per
https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2007-2008/0984.html

the section begins with a series of lists of excluded atoms,
and those lists are then assigned to other atoms by a list
of offsets within that first list... Since it is probably
very unclear, let me use an example:

if there are 8 atoms in total and I want to exclude
{1, 2} from 3 and {1, 4, 5} from 6,
here is the NNB section I need:

5 !NNB
1 2 1 4 5 - column indices
0 0 2 2 2 5 5 5 - row offsets

Basically CSR

*/

// adds explicit nonbond exclusions in psf to kindMap
// pre: stream is before !NNB   post: stream is in groups (NGRP) section just
// after the first appearance of the last donor
//
int ReadPSFExplicitNonbondExclusions(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<unsigned int, std::string>> &firstAtom,
    const uint nNonbondExclusions) {
  Improper imp(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nNonbondExclusions; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &imp.a0, &imp.a1, &imp.a2, &imp.a3);
    if (dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of NNB's in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all NNB's in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this impropers
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign impropers
      if (imp.a0 >= molBegin && imp.a0 < molEnd) {
        imp.a0 -= molBegin;
        imp.a1 -= molBegin;
        imp.a2 -= molBegin;
        imp.a3 -= molBegin;
        // some xplor PSF files have duplicate impropers, we need to ignore
        // these
        if (std::find(currentMol.impropers.begin(), currentMol.impropers.end(),
                      imp) == currentMol.impropers.end()) {
          currentMol.impropers.push_back(imp);
        }
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  return 0;
}

// adds groups in psf to kindMap
// pre: stream is before !NGRP   post: stream is in cross-terms (NCRTERM)
// section just after the first appearance of the last donor
//
int ReadPSFGroups(FILE *psf, MolMap &kindMap,
                  std::vector<std::pair<unsigned int, std::string>> &firstAtom,
                  const uint nGroups) {
  Improper imp(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nGroups; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &imp.a0, &imp.a1, &imp.a2, &imp.a3);
    if (dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of groups in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all groups in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this impropers
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign impropers
      if (imp.a0 >= molBegin && imp.a0 < molEnd) {
        imp.a0 -= molBegin;
        imp.a1 -= molBegin;
        imp.a2 -= molBegin;
        imp.a3 -= molBegin;
        // some xplor PSF files have duplicate impropers, we need to ignore
        // these
        if (std::find(currentMol.impropers.begin(), currentMol.impropers.end(),
                      imp) == currentMol.impropers.end()) {
          currentMol.impropers.push_back(imp);
        }
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  return 0;
}

// adds cross terms in psf to kindMap
// pre: stream is before !NCRTERM   post: stream is at end of psf
// two quadruples of atoms per line:
int ReadPSFCrossTerms(
    FILE *psf, MolMap &kindMap,
    std::vector<std::pair<unsigned int, std::string>> &firstAtom,
    const uint nCrossTerms) {
  Improper imp(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nCrossTerms; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &imp.a0, &imp.a1, &imp.a2, &imp.a3);
    if (dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of cross terms in PSF file ");
      return errors::READ_ERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all cross terms in PSF file ");
      return errors::READ_ERROR;
    }

    // loop to find the molecule kind with this impropers
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind &currentMol = kindMap[firstAtom[i].second];
      // index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      // index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      // assign impropers
      if (imp.a0 >= molBegin && imp.a0 < molEnd) {
        imp.a0 -= molBegin;
        imp.a1 -= molBegin;
        imp.a2 -= molBegin;
        imp.a3 -= molBegin;
        // some xplor PSF files have duplicate impropers, we need to ignore
        // these
        if (std::find(currentMol.impropers.begin(), currentMol.impropers.end(),
                      imp) == currentMol.impropers.end()) {
          currentMol.impropers.push_back(imp);
        }
        // once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  return 0;
}

} // namespace
