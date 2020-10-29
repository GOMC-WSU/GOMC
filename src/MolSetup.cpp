/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "EnsemblePreprocessor.h"
#include "MolSetup.h"
#include "StrLib.h"
#include "ConfigSetup.h"    //For definition of restart
#include "PDBSetup.h"       //For mol names->kinds
#include "FFSetup.h"        //For geometry kinds
#include "BasicTypes.h"
#include "GeomLib.h"

#include <cstdio>
#include <utility>      //for swap (most modern compilers)
#include <algorithm>      //for swap pre-c++11 compilers
#include <cstring>          //strstr

#include <iostream>
#include <iomanip>
#include <sstream>      // std::stringstream



using namespace mol_setup;

bool Dihedral::operator == (const Dihedral& o) const
{
  bool same = false;
  if(a0 == o.a0 && a1 == o.a1 && a2 == o.a2 && a3 == o.a3)
    same = true;

  if(a0 == o.a3 && a1 == o.a2 && a2 == o.a1 && a3 == o.a0)
    same = true;

  return same;
}

bool Dihedral::operator != (const Dihedral& other) const
{
  return !(*this == other);
}

namespace
{
//return if we fail to read anything
const int READERROR = -1;

//Assigns numerical mol kind indices to all molKinds
void AssignMolKinds(MolKind& kind, const pdb_setup::Atoms& pdbData, const std::string& name);
void AssignAtomKinds(MolKind& kind, const FFSetup& ffData);
void AssignBondKinds(MolKind& kind, const FFSetup& ffData);
void AssignAngleKinds(MolKind& kind, const FFSetup& ffData);
void AssignDihKinds(MolKind& kind, const FFSetup& ffData);

void BriefBondKinds(MolKind& kind, const FFSetup& ffData);
void BriefAngleKinds(MolKind& kind, const FFSetup& ffData);
void BriefDihKinds(MolKind& kind, const FFSetup& ffData);

//Builds kindMap from PSF file (does not include coordinates) kindMap
// should be empty returns number of atoms in the file, or READERROR if
// the read failed somehow
int ReadPSF(const char* psfFilename, MolMap& kindMap, SizeMap& sizeMap, pdb_setup::Atoms& pdbData, MolMap * kindMapFromBox1 = NULL, SizeMap * sizeMapFromBox1 = NULL);

//adds atoms and molecule data in psf to kindMap
//pre: stream is at !NATOMS   post: stream is at end of atom section
int ReadPSFAtoms(FILE* psf,
                 MolMap& kindMap, uint nAtoms);
//adds bonds in psf to kindMap
//pre: stream is before !BONDS   post: stream is in bond section just after
//the first appearance of the last molecule
int ReadPSFBonds(FILE* psf, MolMap& kindMap,
                 std::vector<std::pair<uint, std::string> >& firstAtom,
                 const uint nbonds);
//adds angles in psf to kindMap
//pre: stream is before !NTHETA   post: stream is in angle section just
//after the first appearance of the last molecule
int ReadPSFAngles(FILE* psf, MolMap& kindMap,
                  std::vector<std::pair<uint, std::string> >& firstAtom,
                  const uint nangles);
//adds dihedrals in psf to kindMap
//pre: stream is before !NPHI   post: stream is in dihedral section just
//after the first appearance of the last molecule
int ReadPSFDihedrals(FILE* psf, MolMap& kindMap,
                     std::vector<std::pair<uint, std::string> >& firstAtom,
                     const uint ndihedrals);

}


//List of dihedrals with atom at one end, atom first
std::vector<Dihedral> mol_setup::AtomEndDihs(const MolKind& molKind, uint atom)
{
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end(); it < end; ++it) {
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

std::vector<Dihedral> mol_setup::DihsOnBond(const MolKind& molKind, uint atom, uint partner)
{
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end(); it < end; ++it) {
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

std::vector<Dihedral> mol_setup::DihsAll(const MolKind& molKind)
{
  std::vector<Dihedral> result;
  typedef std::vector<Dihedral>::const_iterator Diter;
  for (Diter it = molKind.dihedrals.begin(), end = molKind.dihedrals.end(); it < end; ++it) {
    result.push_back(*it);
  }
  return result;
}

//List of angles with atom at one end, atom first
std::vector<Angle> mol_setup::AtomEndAngles(const MolKind& molKind, uint atom)
{
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end; ++it) {
    if (it->a0 == atom || it->a2 == atom) {
      result.push_back(*it);
    }
    if (it->a2 == atom) {
      std::swap(result.back().a0, result.back().a2);
    }
  }
  return result;
}

std::vector<Angle> mol_setup::AtomMidAngles(const MolKind& molKind, uint atom)
{
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end; ++it) {
    if (it->a1 == atom) {
      result.push_back(*it);
    }
  }
  return result;
}

//List of angles with atom at one end, and mid in middle, atom first
std::vector<Angle> mol_setup::AtomMidEndAngles(const MolKind& molKind, uint mid, uint atom)
{
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end; ++it) {
    if ((it->a0 == atom || it->a2 == atom) && (it->a1 == mid)) {
      result.push_back(*it);

      if (it->a2 == atom) {
        std::swap(result.back().a0, result.back().a2);
      }
    }
  }
  return result;
}

std::vector<Angle> mol_setup::AngsAll(const MolKind& molKind)
{
  std::vector<Angle> result;
  typedef std::vector<Angle>::const_iterator Aiter;
  for (Aiter it = molKind.angles.begin(), end = molKind.angles.end(); it < end; ++it) {
    result.push_back(*it);
  }
  return result;
}

//List of bonds with atom at one end, atom first
std::vector<Bond> mol_setup::AtomBonds(const MolKind& molKind, uint atom)
{
  std::vector<Bond> result;
  typedef std::vector<Bond>::const_iterator Biter;
  for (Biter it = molKind.bonds.begin(), end = molKind.bonds.end(); it < end; ++it) {
    if (it->a0 == atom || it->a1 == atom) {
      result.push_back(*it);
    }
    if (it->a1 == atom) {
      std::swap(result.back().a0, result.back().a1);
    }
  }
  return result;
}

std::vector<Bond> mol_setup::BondsAll(const MolKind& molKind)
{
  std::vector<Bond> result;
  typedef std::vector<Bond>::const_iterator Biter;
  for (Biter it = molKind.bonds.begin(), end = molKind.bonds.end(); it < end; ++it) {
    result.push_back(*it);
  }
  return result;
}

int mol_setup::ReadCombinePSF(MolMap& kindMap,
                              SizeMap& sizeMap,
                              std::string const*const psfFilename,
                              const int numFiles, pdb_setup::Atoms& pdbAtoms)
{
  int errorcode = ReadPSF(psfFilename[0].c_str(), kindMap, sizeMap, pdbAtoms);
  if (errorcode < 0)
    return errorcode;
  MolMap map2;
  SizeMap sizeMap2;
  for (int i = 1; i < numFiles; ++i) {
    map2.clear();
    errorcode = ReadPSF(psfFilename[i].c_str(), map2, sizeMap2, pdbAtoms, &kindMap, &sizeMap);

    if (errorcode < 0)
      return errorcode;
    //kindMap.insert(map2.begin(), map2.end());
  }

  PrintMolMapVerbose(kindMap);
  //PrintMolMapBrief(kindMap);

  return 0;
}
int MolSetup::Init(const config_setup::RestartSettings& restart,
                   const std::string* psfFilename, pdb_setup::Atoms& pdbAtoms)
{
  kindMap.clear();
  sizeMap.clear();

  int numFiles;
  if(restart.enable)
    numFiles = 1;
  else
    numFiles = BOX_TOTAL;

  return ReadCombinePSF(kindMap, sizeMap, psfFilename, numFiles, pdbAtoms);

}


void MolSetup::AssignKinds(const pdb_setup::Atoms& pdbAtoms, const FFSetup& ffData)
{
  typedef MolMap::iterator MapIt;
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    AssignMolKinds(it->second, pdbAtoms, it->first);
    AssignAtomKinds(it->second, ffData);
    AssignBondKinds(it->second, ffData);
    AssignAngleKinds(it->second, ffData);
    AssignDihKinds(it->second, ffData);
  }

  //Print bonded Information
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
  printf("%s %33s %4s %16s \n", "Atom Types", "Kchi(K)", "n",
         "delta(degree)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefDihKinds(it->second, ffData);
  }
  std::cout << std::endl;
}

int read_atoms(FILE *psf, unsigned int nAtoms, std::vector<mol_setup::Atom> & allAtoms)
{
  char input[512];
  unsigned int atomID = 0;
  unsigned int molID;
  char segment[11], moleculeName[11], atomName[11], atomType[11];
  double charge, mass;

  while (atomID < nAtoms) {
    char* check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Could not find all atoms in PSF file ");
      return READERROR;
    }
    //skip comment/blank lines
    if (input[0] == '!' || str::AllWS(input))
      continue;
    //parse line
    sscanf(input, " %u %s %u %s %s %s %lf %lf ",
           &atomID, segment, &molID,
           moleculeName, atomName, atomType, &charge, &mass);
    allAtoms.push_back(mol_setup::Atom(atomName, moleculeName, molID, segment, atomType, charge, mass));
  }
}

typedef std::vector<uint>::const_iterator candidateIterator;

int createMapAndModifyPDBAtomDataStructure( const BondAdjacencyList & bondAdjList,
                                            const std::vector< std::vector<uint> > & moleculeXAtomIDY, 
                                            std::vector<mol_setup::Atom> & allAtoms,
                                            mol_setup::MolMap & kindMap,
                                            mol_setup::SizeMap & sizeMap,
                                            pdb_setup::Atoms& pdbAtoms,
                                            mol_setup::MolMap * kindMapFromBox1,
                                            mol_setup::SizeMap * sizeMapFromBox1){

  /* A size -> moleculeKind map for quick evaluation of new molecules based on molMap entries
    of a given size exisitng or not */ 
  uint startIdxResBoxOffset;
  uint resKindIndex;
  if (pdbAtoms.lastAtomIndexInBox0 == 0){
  startIdxResBoxOffset = 0;
  resKindIndex = 0;
  pdbAtoms.startIdxRes.clear();
  pdbAtoms.resKinds.clear();
  pdbAtoms.resKindNames.clear();
  pdbAtoms.resNames.clear();
  } else {
    startIdxResBoxOffset = pdbAtoms.lastAtomIndexInBox0 + 1;
    resKindIndex = pdbAtoms.lastResKindIndex;
  }


  /* Iterate through N connected components */
  int stringSuffix = 1;
  for (std::vector< std::vector<uint> >::const_iterator it = moleculeXAtomIDY.cbegin();
        it != moleculeXAtomIDY.cend(); it++){

    std::string fragName;
    bool multiResidue = false;
    bool newSize = false;
    bool newMapEntry = true;
    bool foundEntryInOldMap = true;

    /* Search by size for existing molecules from Box 1 if it exists*/
    if (sizeMapFromBox1 != NULL && kindMapFromBox1 != NULL){
      SizeMap::iterator sizeIt = sizeMapFromBox1->find(it->size());

      /* Found no matching molecules from Box 1 by size */
      if (sizeIt == sizeMapFromBox1->end()) {
        /* For record keeping later on */
        foundEntryInOldMap = false;
      /* Found molecules with the same size, now evaluate for atom equality */
      } else {
        /* Iterate through all the size consistent map entries */
        for (std::vector<std::__cxx11::string>::const_iterator sizeConsistentEntries = sizeIt->second.cbegin();
          sizeConsistentEntries != sizeIt->second.cend(); sizeConsistentEntries++){
          /* Iterate atom by atom of a given size consistent map entries with the candidate molecule*/
          typedef std::vector<mol_setup::Atom>::const_iterator atomIterator;
          std::pair<atomIterator, candidateIterator> itPair((*kindMapFromBox1)[*sizeConsistentEntries].atoms.cbegin(), it->cbegin());
          for (; itPair.second != it->cend(); ++itPair.first, ++itPair.second){
            if (itPair.first->name == allAtoms[*(itPair.second)].name){
              continue;
            } else {
              break;
            }
          }
          // Found a match in our our molMap.
          if (itPair.second == it->cend()) {

            /* Get the map key */
            fragName = *sizeConsistentEntries;
            /* Boilerplate PDB Data modifications for matches */
            pdbAtoms.startIdxRes.push_back(startIdxResBoxOffset + it->front());
            pdbAtoms.resKinds.push_back((*kindMapFromBox1)[fragName].kindIndex);
            pdbAtoms.resNames.push_back(fragName);
            newMapEntry = false;
            /* Boilerplate PDB Data modifications for matches */

            /* Search current KindMap for this entry. 
              We won't have a value for fragName unless we matched
              in the old molMap, since it may be a PROTX name.
              That's why this all has to be in this conditional.
            */ 
            MolMap::const_iterator kindIt = kindMap.find(fragName);

            /* If we don't have it in our new Map, then we need to add the first 
              entry in our new Map to quench the post processessing requirements. */
            if (kindIt == kindMap.cend()){
              if((*kindMapFromBox1)[fragName].isMultiResidue){
                kindMap[fragName] = MolKind();
                kindMap[fragName].isMultiResidue = true;
                for (std::vector<uint>::const_iterator connectedComponentIt = it->cbegin();
                connectedComponentIt != it->cend(); connectedComponentIt++){
                  kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
                  kindMap[fragName].intraMoleculeResIDs.push_back(allAtoms[*connectedComponentIt].residueID);
                }
                /* Normalize resID to intramolecule indices */
                uint firstResID = kindMap[fragName].intraMoleculeResIDs.front();
                for (auto& x : kindMap[fragName].intraMoleculeResIDs)
                  x -= firstResID;
                /* Normalize resID to intramolecule indices */
              } else {
                kindMap[fragName] = MolKind();
                kindMap[fragName].isMultiResidue = false;
                for (std::vector<uint>::const_iterator connectedComponentIt = it->cbegin();
                connectedComponentIt != it->cend(); connectedComponentIt++){
                  kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
                }
              }
              kindMap[fragName].firstAtomID = it->front() + 1;
              kindMap[fragName].firstMolID = allAtoms[it->front()].residueID;
              kindMap[fragName].kindIndex = (*kindMapFromBox1)[fragName].kindIndex;
              MolSetup::copyBondInfoIntoMapEntry(bondAdjList, kindMap, fragName);

              /* Finally, search new sizeMap for existing entry of same size as old molecule.
                 This handles chance that we have two equal sized, but different, molecules in
                 the two boxes.  Since the map entry will already exist, it wouldn't be a newSize */

              SizeMap::iterator sizeIt = sizeMap.find(it->size());
              /* New Size */
              if (sizeIt == sizeMap.end()) {
                sizeMap[it->size()] = std::vector<std::__cxx11::string>{fragName};
              } else {
                sizeMap[it->size()].push_back(fragName);
              }
            }

          } else {
            foundEntryInOldMap = false;
          }
        }
      }
    } else {
      foundEntryInOldMap = false;
    }

    /* Follow procedure from earlier versions of single box ensemble map entries. Very similar to above,
       except we may need to generate fragName in the case of a protein. */
    if(!foundEntryInOldMap){

      /* Search by size for existing molecules  for quick evaluation of new molecules based on molMap entries
      of a given size exisitng or not */ 
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
        for (std::vector<std::__cxx11::string>::const_iterator sizeConsistentEntries = sizeIt->second.cbegin();
          sizeConsistentEntries != sizeIt->second.cend(); sizeConsistentEntries++){
          /* Iterate atom by atom of a given size consistent map entries with the candidate molecule*/
          typedef std::vector<mol_setup::Atom>::const_iterator atomIterator;
          std::pair<atomIterator, candidateIterator> itPair(kindMap[*sizeConsistentEntries].atoms.cbegin(), it->cbegin());
          for (; itPair.second != it->cend(); ++itPair.first, ++itPair.second){
            if (itPair.first->name == allAtoms[*(itPair.second)].name){
              continue;
            } else {
              break;
            }
          }
          // Found a match
          if (itPair.second == it->cend()) {
            // Modify PDBData
            pdbAtoms.startIdxRes.push_back(startIdxResBoxOffset + it->front());
            pdbAtoms.resKinds.push_back(kindMap[*sizeConsistentEntries].kindIndex);
            newMapEntry = false;
          }
        }
      }

      if (newMapEntry){
        /* Determine if this connected component is a standard or multiResidue molecule */
        for (candidateIterator connectedComponentIt = it->cbegin();
          connectedComponentIt != it->cend(); connectedComponentIt++){
          if (allAtoms[*connectedComponentIt].residueID == allAtoms[it->front()].residueID){
            continue;
          } else {
            multiResidue = true;
            break;
          }
        }
      }

      if (newMapEntry){
        if(multiResidue){  
          std::stringstream ss;
          /* Length of Suffix */
          for (int i = 0; i < (stringSuffix + 26 - 1) / 26; i++){
            int intermediate = stringSuffix - i * 26;
            char charSuffix = 'A';
            /* Increment Char A until reach suffix or 27 which will be Z. */
            for (int j = 1; j < std::min(intermediate, 27); j++){
              charSuffix++;
            }
            ss << charSuffix;
          }
          fragName = "PROT" + ss.str();
          printf("\n%-40s \n", "Warning: A molecule containing > 1 residue is detected.");
          printf("The simulation will name it %s.\n", fragName.c_str());
          printf("See the chart at the end of the output log describing this entry.\n");

          kindMap[fragName] = MolKind();
          kindMap[fragName].isMultiResidue = true;
          for (std::vector<uint>::const_iterator connectedComponentIt = it->cbegin();
          connectedComponentIt != it->cend(); connectedComponentIt++){
            kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
            kindMap[fragName].intraMoleculeResIDs.push_back(allAtoms[*connectedComponentIt].residueID);
          }
          /* Normalize resID to intramolecule indices */
          uint firstResID = kindMap[fragName].intraMoleculeResIDs.front();
          for (auto& x : kindMap[fragName].intraMoleculeResIDs)
            x -= firstResID;
          /* Normalize resID to intramolecule indices */
        } else {
          fragName = allAtoms[it->front()].residue;
          kindMap[allAtoms[it->front()].residue] = MolKind();
          kindMap[fragName].isMultiResidue = false;
          for (std::vector<uint>::const_iterator connectedComponentIt = it->cbegin();
          connectedComponentIt != it->cend(); connectedComponentIt++){
            kindMap[fragName].atoms.push_back(allAtoms[*connectedComponentIt]);
          }
        }
        kindMap[fragName].firstAtomID = it->front() + 1;
        kindMap[fragName].firstMolID = allAtoms[it->front()].residueID;
        kindMap[fragName].kindIndex = resKindIndex;
        //pdbAtoms.startIdxRes.push_back(kindMap[fragName].firstAtomID - 1);
        pdbAtoms.startIdxRes.push_back(startIdxResBoxOffset + kindMap[fragName].firstAtomID - 1);
        pdbAtoms.resKinds.push_back(kindMap[fragName].kindIndex);
        pdbAtoms.resKindNames.push_back(fragName);
        pdbAtoms.resNames.push_back(fragName);
        MolSetup::copyBondInfoIntoMapEntry(bondAdjList, kindMap, fragName);
        resKindIndex++;
        if (newSize){
          sizeMap[it->size()] = std::vector<std::__cxx11::string>{fragName};
        } else {
          sizeMap[it->size()].push_back(fragName);
        }
      }
    }
  }

  pdbAtoms.lastAtomIndexInBox0 = (moleculeXAtomIDY.back()).back();
  pdbAtoms.lastResKindIndex = resKindIndex;
}

typedef std::map<std::__cxx11::string, mol_setup::MolKind> MolMap;
void MolSetup::copyBondInfoIntoMapEntry(const BondAdjacencyList & bondAdjList, mol_setup::MolMap & kindMap, std::string fragName){

    unsigned int molBegin = kindMap[fragName].firstAtomID - 1;
    //index AFTER last atom in molecule
    unsigned int molEnd = molBegin + kindMap[fragName].atoms.size();
    //assign the bond
    for (uint i = molBegin; i < molEnd; i++){
      adjNode* ptr = bondAdjList.head[i];
      while (ptr != nullptr) {
        if (i < ptr->val) {
          kindMap[fragName].bonds.push_back(Bond(i-molBegin, ptr->val-molBegin));
        }    
        ptr = ptr->next;
      }
    }
}

namespace
{

void AssignMolKinds(MolKind& kind, const pdb_setup::Atoms& pdbData, const std::string& name)
{
  /* Bug in old code, should subtract beginning not end */
  uint index = std::find(pdbData.resKindNames.begin(),
                         pdbData.resKindNames.end(), name) - pdbData.resKindNames.begin();
  kind.kindIndex = index;
}

void AssignAtomKinds(MolKind& kind, const FFSetup& ffData)
{
  for (uint i = 0; i < kind.atoms.size(); ++i) {
    int thisKind = ffData.mie.Find(&kind.atoms[i].type, ffData.mie.name);
    if (thisKind < 0) {
      fprintf(stderr, "ERROR: Atom Type %s not specified in nonbonded section of parameter file.\n", kind.atoms[i].type.c_str());
      exit(EXIT_FAILURE);
    }
    kind.atoms[i].kind = thisKind;
  }
}


void AssignBondKinds(MolKind& kind, const FFSetup& ffData)
{
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

void BriefBondKinds(MolKind& kind, const FFSetup& ffData)
{
  const uint ATOMS_PER = 2;
  std::string elementNames[ATOMS_PER];
  std::vector<std::string> printed;

  if(kind.bonds.size() == 0)
    return;

  for(uint i = 0; i < kind.bonds.size(); ++i) {
    uint search = kind.bonds[i].kind;
    std::string bondName;

    elementNames[0] = kind.atoms[kind.bonds[i].a0].type;
    elementNames[1] = kind.atoms[kind.bonds[i].a1].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      bondName.append(elementNames[m]).append("\t");
    }

    if(find(printed.begin(), printed.end(), bondName) == printed.end()) {
      printf("%s", bondName.c_str());
      if(ffData.bond.GetKb(search) > 99999999)
        printf("%28s %16.4f \n", "FIX", ffData.bond.Getb0(search));
      else
        printf("%28.4f %16.4f \n", ffData.bond.GetKb(search),
               ffData.bond.Getb0(search));

      printed.push_back(bondName);
    }
  }
  std::cout << std::endl;
}

void AssignAngleKinds(MolKind& kind, const FFSetup& ffData)
{
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

void BriefAngleKinds(MolKind& kind, const FFSetup& ffData)
{
  const uint ATOMS_PER = 3;
  std::string elementNames[ATOMS_PER];
  std::vector<std::string> printed;
  double coef = 180.00 / M_PI;

  if(kind.angles.size() == 0)
    return;

  for(uint i = 0; i < kind.angles.size(); ++i) {
    std::string angleName;
    uint search = kind.angles[i].kind;
    elementNames[0] = kind.atoms[kind.angles[i].a0].type;
    elementNames[1] = kind.atoms[kind.angles[i].a1].type;
    elementNames[2] = kind.atoms[kind.angles[i].a2].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      angleName.append(elementNames[m]).append("\t");
    }

    if(find(printed.begin(), printed.end(), angleName) == printed.end()) {
      printf("%s", angleName.c_str());
      if(ffData.angle.GetKtheta(search) > 99999999)
        printf("%20s %16.4f \n", "FIX", ffData.angle.Gettheta0(search) *coef);
      else
        printf("%20.4f %16.4f \n", ffData.angle.GetKtheta(search),
               ffData.angle.Gettheta0(search) * coef);

      printed.push_back(angleName);
    }
  }
  std::cout << std::endl;
}

void AssignDihKinds(MolKind& kind, const FFSetup& ffData)
{
  const uint ATOMS_PER = 4;
  std::string elementNames[ATOMS_PER];

  int search = 0;
  for(uint i = 0; i < kind.dihedrals.size(); ++i) {
    elementNames[0] = kind.atoms[kind.dihedrals[i].a0].type;
    elementNames[1] = kind.atoms[kind.dihedrals[i].a1].type;
    elementNames[2] = kind.atoms[kind.dihedrals[i].a2].type;
    elementNames[3] = kind.atoms[kind.dihedrals[i].a3].type;
    search = ffData.dih.Find(elementNames, ffData.dih.name);
    if(search < 0) {
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

void BriefDihKinds(MolKind& kind, const FFSetup& ffData)
{
  const uint ATOMS_PER = 4;
  std::string elementNames[ATOMS_PER];
  double coef = 180.00 / M_PI;
  std::vector<std::string> printed;

  if(kind.dihedrals.size() == 0)
    return;

  for(uint i = 0; i < kind.dihedrals.size(); ++i) {
    std::string dName = ffData.dih.name[kind.dihedrals[i].kind];
    std::string dihedralName;
    uint dihsize = ffData.dih.GetSizeDih(dName);

    elementNames[0] = kind.atoms[kind.dihedrals[i].a0].type;
    elementNames[1] = kind.atoms[kind.dihedrals[i].a1].type;
    elementNames[2] = kind.atoms[kind.dihedrals[i].a2].type;
    elementNames[3] = kind.atoms[kind.dihedrals[i].a3].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      dihedralName.append(elementNames[m]).append("\t");
    }

    if(find(printed.begin(), printed.end(), dihedralName) == printed.end()) {
      for(uint j = 0; j < dihsize; j++) {
        printf("%s", dihedralName.c_str());
        printf("%12.4f %4d %11.4f \n", ffData.dih.GetKchi(dName, j),
               ffData.dih.Getn(dName, j),
               ffData.dih.Getdelta(dName, j) * coef);
      }
      printed.push_back(dihedralName);
    }
  }
}

}



void mol_setup::PrintMolMapVerbose(const MolMap& kindMap)
{
  std::cout << "\nMolecules in PSF:\n";
  MolMap::const_iterator it = kindMap.begin();
  while (it != kindMap.end()) {
    std::cout << "Molecule Kind: " << it->first << std::endl;
    std::cout << "Idx\tname\ttype\tcharge\tmass\n";
    for (uint i = 0; i < it->second.atoms.size(); i++) {
      std::cout << i << "\t" << it->second.atoms[i].name << '\t' <<
                it->second.atoms[i].type << '\t' << std::setprecision(4) <<
                it->second.atoms[i].charge << '\t' << std::setprecision(4) <<
                it->second.atoms[i].mass << std::endl;
    }
    std::cout << "\nBonds:";
    for (uint i = 0; i < it->second.bonds.size(); i++) {
      if (i % 10 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.bonds[i].a0 << ' ' << it->second.bonds[i].a1 << ']' << ' ';
    }
    std::cout << std::endl << "\nAngles:";
    for (uint i = 0; i < it->second.angles.size(); i++) {
      if (i % 7 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.angles[i].a0 << ' '
                << it->second.angles[i].a1 << ' '
                << it->second.angles[i].a2 << ']' << ' ';
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
    ++it;
    std::cout << std::endl << std::endl;
  }
}

void mol_setup::PrintMolMapBrief(const MolMap& kindMap)
{
  std::cout << "Molecules in PSF:\n";
  std::cout << "Name\t#Atom\t#Bond\t#Ang\t#Dih\t\n";
  MolMap::const_iterator it = kindMap.begin();
  while (it != kindMap.end()) {
    std::cout << it->first << '\t' << it->second.atoms.size() << '\t' <<
              it->second.bonds.size() << '\t' <<
              it->second.angles.size() << '\t' <<
              it->second.dihedrals.size() << '\t' << std::endl;
    ++it;
  }

}


namespace
{
//Initializes system from PSF file (does not include coordinates)
//returns number of atoms in the file, or READERROR if the read failed somehow
int ReadPSF(const char* psfFilename, MolMap& kindMap, SizeMap & sizeMap, pdb_setup::Atoms& pdbData, MolMap * kindMapFromBox1, SizeMap * sizeMapFromBox1)
{
  FILE* psf = fopen(psfFilename, "r");
  char* check;		//return value of fgets
  int count;		//for number of bonds/angles/dihs
  if (psf == NULL) {
    fprintf(stderr, "ERROR: Failed to open PSF file %s for molecule data.\nExiting...\n", psfFilename);
    return READERROR;
  }
  char input[512];
  unsigned int nAtoms;
  //find atom header+count
  do {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read atoms from PSF file %s",
              psfFilename);
      fclose(psf);
      return READERROR;
    }
  } while (strstr(input, "!NATOM") == NULL);
  sscanf(input, " %u", &nAtoms);

  /* GJS - This is a flat vector of atom objects.  Instead of building kindMap entries as we parse
  the PSF file, we will first read the entire atom section, building an N length array of atoms,
  N is number of atoms.  Since, we process the PSF file as a flat file handle, we cannot process
  bonds before atoms without physically generating a new PSFFile reversing the order of ATOMS <-> BONDS.
  Hence, the necessity to build this vector before knowing how the atoms are connected. */
  std::vector<mol_setup::Atom> allAtoms;
  read_atoms(psf, nAtoms, allAtoms);
  //build list of start particles for each type, so we can find it and skip
  //everything else
  //make sure molecule has bonds, appears before !NBOND
    //find bond header+count
  while (strstr(input, "!NBOND") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read bonds from PSF file %s",
              psfFilename);
      fclose(psf);
      return  READERROR;
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
  std::vector< std::vector<uint> > moleculeXAtomIDY;
  /* A standard adjacency list with N nodes, where N is number of atoms.
   This is an undirected graph, where edges between nodes
   represent bonds between atoms.  It is generated by DFS, checking if a node
   has been visited before.
  */
  BondAdjacencyList bondAdjList(psf, nAtoms, count, moleculeXAtomIDY);
  /* Molecular equality is determined by a series of evaluations, with early termination - 
    1) Length of candidate versus all entries in the map
    2) Atom by atom check for equality between candidate and all entries.

    If the candidate matches an exisiting map entry, we push the first resID onto startIDxRes, along
    with other necessary information.

    If the candidate is determined to be novel, we first determine its status as standard or multiresidue.
    
    If multiresidue a dummy string is generated and atoms added to the entry along with intraMolecularResIDs
    normalized to 0 .. n-1, n is number of residues in the multiresidue molecule.

    Otherwise, standard procedure for creating a map entry is followed.

    The bond information contained in the Adjacency list is assigned to map entries.

    Finally, entries in startIDxRes are consolidated redefine the start and end of molecule,
    as far as the pdb data is concerned. 
  */
  createMapAndModifyPDBAtomDataStructure( bondAdjList, 
                                          moleculeXAtomIDY, 
                                          allAtoms, 
                                          kindMap, 
                                          sizeMap,
                                          pdbData, 
                                          kindMapFromBox1, 
                                          sizeMapFromBox1);

  std::vector<std::pair<unsigned int, std::string> > firstAtomLookup;
  for (MolMap::iterator it = kindMap.begin(); it != kindMap.end(); ++it) {
    firstAtomLookup.push_back(std::make_pair(it->second.firstAtomID, it->first));
  }
  std::sort(firstAtomLookup.begin(), firstAtomLookup.end());
  //find bond header+count
  //make sure molecule has bonds, appears before !NBOND
  //find angle header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NTHETA") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read angles from PSF file %s",
              psfFilename);
      fclose(psf);
      return READERROR;
    }
  }
  //make sure molecule has angles, count appears before !NTHETA
  count = atoi(input);
  if (ReadPSFAngles(psf, kindMap, firstAtomLookup, count) == READERROR) {
    fclose(psf);
    return READERROR;
  }
  //find dihedrals header+count
  fseek(psf, 0, SEEK_SET);
  while (strstr(input, "!NPHI") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read dihedrals from PSF file %s",
              psfFilename);
      fclose(psf);
      return READERROR;
    }
  }
  //make sure molecule has dihs, count appears before !NPHI
  count = atoi(input);

  if (ReadPSFDihedrals(psf, kindMap, firstAtomLookup, count) == READERROR) {
    fclose(psf);
    return READERROR;
  }

  fclose(psf);

  return nAtoms;
}

//adds atoms and molecule data in psf to kindMap
//pre: stream is at !NATOMS   post: stream is at end of atom section
int ReadPSFAtoms(FILE* psf, MolMap& kindMap, unsigned int nAtoms)
{
  char input[512];
  unsigned int atomID = 0;
  unsigned int molID;
  char segment[11], moleculeName[11], atomName[11], atomType[11];
  double charge, mass;

  while (atomID < nAtoms) {
    char* check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Could not find all atoms in PSF file ");
      return READERROR;
    }
    //skip comment/blank lines
    if (input[0] == '!' || str::AllWS(input))
      continue;
    //parse line
    sscanf(input, " %u %s %u %s %s %s %lf %lf ",
           &atomID, segment, &molID,
           moleculeName, atomName, atomType, &charge, &mass);
    MolMap::iterator it = kindMap.find(moleculeName);
    //found new molecule kind...
    if (it == kindMap.end()) {
      it = kindMap.insert(std::make_pair(std::string(moleculeName), MolKind())).first;
      it->second.firstAtomID = atomID;
      it->second.firstMolID = molID;
      it->second.atoms.push_back(Atom(atomName, moleculeName, molID, segment, atomType, charge, mass));
    }
    //still building a molecule...
    else if (it->second.incomplete) {
      if (molID != it->second.firstMolID)
        it->second.incomplete = false;
      else
        it->second.atoms.push_back(Atom(atomName, moleculeName, molID, segment, atomType, charge, mass));
    }
  }
  //Fix for one molecule fringe case.
  if (molID == 1) {
    MolMap::iterator it = kindMap.find(moleculeName);
    it->second.incomplete = false;
  }
  return 0;
}

//adds bonds in psf to kindMap
//pre: stream is before !BONDS   post: stream is in bond section just after
//the first appearance of the last molecule
int ReadPSFBonds(FILE* psf, MolMap& kindMap,
                 std::vector<std::pair<unsigned int, std::string> >& firstAtom,
                 const uint nbonds)
{
  unsigned int atom0, atom1;
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nbonds; n++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if(dummy != 2) {
      fprintf(stderr, "ERROR: Incorrect Number of bonds in PSF file ");
      return READERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all bonds in PSF file ");
      return READERROR;
    }

    //loop to find the molecule kind with this bond
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind& currentMol = kindMap[firstAtom[i].second];
      //index of first atom in molecule
      unsigned int molBegin = firstAtom[i].first;
      //index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      //assign the bond
      if (atom0 >= molBegin && atom0 < molEnd) {
        currentMol.bonds.push_back(Bond(atom0 - molBegin, atom1 - molBegin));
        //once we found the molecule kind, break from the loop
        defined[i] = true;
        break;
      }
    }
  }
  //Check if we defined all bonds
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    if(currentMol.atoms.size() > 1 && !defined[i]) {
      std::cout << "Warning: Bond is missing for " << firstAtom[i].second
                << " !\n";
    }
  }
  return 0;
}

//adds angles in psf to kindMap
//pre: stream is before !NTHETA   post: stream is in angle section just after
//the first appearance of the last molecule
int ReadPSFAngles(FILE* psf, MolMap& kindMap,
                  std::vector<std::pair<unsigned int, std::string> >& firstAtom,
                  const uint nangles)
{
  unsigned int atom0, atom1, atom2;
  int dummy;
  //std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < nangles; n++) {
    dummy = fscanf(psf, "%u %u %u", &atom0, &atom1, &atom2);
    if(dummy != 3) {
      fprintf(stderr, "ERROR: Incorrect Number of angles in PSF file ");
      return READERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all angles in PSF file ");
      return READERROR;
    }

    //loop to find the molecule kind with this angles
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind& currentMol = kindMap[firstAtom[i].second];
      //index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      //index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      //assign the angle
      if (atom0 >= molBegin && atom0 < molEnd) {
        currentMol.angles.push_back(Angle(atom0 - molBegin, atom1 - molBegin,
                                          atom2 - molBegin));
        //once we found the molecule kind, break from the loop
        currentMol.anglesDefined = true;
        break;
      }
    }
  }
  //Check if we defined all angles
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    if(currentMol.atoms.size() > 2 && !currentMol.anglesDefined) {
      std::cout << "Warning: Angle is missing for " << firstAtom[i].second
                << " !\n";
    }
  }
  return 0;
}



bool ContainsDihedral(const std::vector<uint>& vec, const int dih[])
{
  for (uint i = 0; i < vec.size(); i += 4) {
    bool match = true;
    for (uint j = 0; j < 4; ++j)
      if (vec[i + j] != dih[j])
        match = false;
    if (match)
      return true;
  }
  return false;
}


//adds dihedrals in psf to kindMap
//pre: stream is before !NPHI   post: stream is in dihedral section just after
//the first appearance of the last molecule
//
int ReadPSFDihedrals(FILE* psf, MolMap& kindMap,
                     std::vector<std::pair<unsigned int, std::string> >& firstAtom, const uint ndihedrals)
{
  Dihedral dih(0, 0, 0, 0);
  int dummy;
  std::vector<bool> defined(firstAtom.size(), false);
  for (uint n = 0; n < ndihedrals; n++) {
    dummy = fscanf(psf, "%u %u %u %u", &dih.a0, &dih.a1, &dih.a2, &dih.a3);
    if(dummy != 4) {
      fprintf(stderr, "ERROR: Incorrect Number of dihedrals in PSF file ");
      return READERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all dihedrals in PSF file ");
      return READERROR;
    }

    //loop to find the molecule kind with this dihedral
    for (unsigned int i = 0; i < firstAtom.size(); ++i) {
      MolKind& currentMol = kindMap[firstAtom[i].second];
      //index of first atom in moleule
      unsigned int molBegin = firstAtom[i].first;
      //index AFTER last atom in molecule
      unsigned int molEnd = molBegin + currentMol.atoms.size();
      //assign dihedral
      if (dih.a0 >= molBegin && dih.a0 < molEnd) {
        dih.a0 -= molBegin;
        dih.a1 -= molBegin;
        dih.a2 -= molBegin;
        dih.a3 -= molBegin;
        //some xplor PSF files have duplicate dihedrals, we need to ignore these
        if (std::find(currentMol.dihedrals.begin(), currentMol.dihedrals.end(),
                      dih) == currentMol.dihedrals.end()) {
          currentMol.dihedrals.push_back(dih);
        }
        //once we found the molecule kind, break from the loop
        currentMol.dihedralsDefined = true;
        break;
      }
    }
  }
  //Check if we defined all dihedrals in the map so far.
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    if(currentMol.atoms.size() > 3 && !currentMol.dihedralsDefined) {
      std::cout << "Warning: Dihedral is missing for " << firstAtom[i].second
                << " !\n";
    }
  }
  return 0;
}

}
