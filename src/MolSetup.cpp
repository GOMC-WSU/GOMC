/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
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



using namespace mol_setup;

bool Dihedral::operator == (const Dihedral& o) const
{
  return (a0 == o.a0 && a1 == o.a1 && a2 == o.a2 && a3 == o.a3);
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
int ReadPSF(const char* psfFilename, MolMap& kindMap);
//adds atoms and molecule data in psf to kindMap
//pre: stream is at !NATOMS   post: stream is at end of atom section
int ReadPSFAtoms(FILE* psf,
                 MolMap& kindMap, uint nAtoms);
//adds bonds in psf to kindMap
//pre: stream is before !BONDS   post: stream is in bond section just after
//the first appearance of the last molecule
int ReadPSFBonds(FILE* psf, MolMap& kindMap,
                 std::vector<std::pair<uint, std::string> >& firstAtom);
//adds angles in psf to kindMap
//pre: stream is before !NTHETA   post: stream is in angle section just
//after the first appearance of the last molecule
int ReadPSFAngles(FILE* psf, MolMap& kindMap,
                  std::vector<std::pair<uint, std::string> >& firstAtom);
//adds dihedrals in psf to kindMap
//pre: stream is before !NPHI   post: stream is in dihedral section just
//after the first appearance of the last molecule
int ReadPSFDihedrals(FILE* psf, MolMap& kindMap,
                     std::vector<std::pair<uint, std::string> >& firstAtom);

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

int mol_setup::ReadCombinePSF(MolMap& kindMap,
                              std::string const*const psfFilename,
                              const int numFiles)
{
  int errorcode = ReadPSF(psfFilename[0].c_str(), kindMap);
  if (errorcode < 0)
    return errorcode;
  MolMap map2;
  for (int i = 1; i < numFiles; ++i) {
    map2.clear();
    errorcode = ReadPSF(psfFilename[i].c_str(), map2);
    if (errorcode < 0)
      return errorcode;
    kindMap.insert(map2.begin(), map2.end());
  }

  PrintMolMapVerbose(kindMap);
  //PrintMolMapBrief(kindMap);

  return 0;
}
int MolSetup::Init(const config_setup::RestartSettings& restart,
                   const std::string* psfFilename)
{
  kindMap.clear();
  int numFiles;
  if(restart.enable)
    numFiles = 1;
  else
    numFiles = BOX_TOTAL;

  return ReadCombinePSF(kindMap, psfFilename, numFiles);
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
  printf("%-19s %15s %20s \n", "Atom Types", "Kb(K)", "b0(A)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefBondKinds(it->second, ffData);
  }

  printf("Angles parameter:\n");
  printf("%-19s %15s %20s \n", "Atom Types", "Ktheta(K)", "theta0(degree)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefAngleKinds(it->second, ffData);
  }

  printf("Dihedrals parameter:\n");
  printf("%-19s %15s %4s %15s \n", "Atom Types", "Kchi(K)", "n",
         "delta(degree)");
  for (MapIt it = kindMap.begin(), end = kindMap.end(); it != end; ++it) {
    BriefDihKinds(it->second, ffData);
  }
  std::cout << std::endl;
}

namespace
{

void AssignMolKinds(MolKind& kind, const pdb_setup::Atoms& pdbData, const std::string& name)
{
  uint index = std::find(pdbData.resKindNames.begin(),
                         pdbData.resKindNames.end(), name) - pdbData.resKindNames.end();
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
    std::string bondName, bondNameReverse;

    elementNames[0] = kind.atoms[kind.bonds[i].a0].type;
    elementNames[1] = kind.atoms[kind.bonds[i].a1].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      bondName.append(elementNames[m]).append("  ");
      bondNameReverse.append(elementNames[ATOMS_PER - m - 1]).append("  ");
    }

    if(find(printed.begin(), printed.end(), bondName) == printed.end()) {
      printf("%-20s", bondName.c_str());
      if(ffData.bond.GetKb(search) > 99999999)
        printf("%15s %20.4f \n", "FIX", ffData.bond.Getb0(search));
      else
        printf("%15.6f %20.4f \n", ffData.bond.GetKb(search),
               ffData.bond.Getb0(search));

      printed.push_back(bondName);
      printed.push_back(bondNameReverse);
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
    std::string angleName, angleNameReverse;
    uint search = kind.angles[i].kind;
    elementNames[0] = kind.atoms[kind.angles[i].a0].type;
    elementNames[1] = kind.atoms[kind.angles[i].a1].type;
    elementNames[2] = kind.atoms[kind.angles[i].a2].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      angleName.append(elementNames[m]).append("  ");
      angleNameReverse.append(elementNames[ATOMS_PER - m - 1]).append("  ");
    }

    if(find(printed.begin(), printed.end(), angleName) == printed.end()) {
      printf("%-20s", angleName.c_str());
      if(ffData.angle.GetKtheta(search) > 99999999)
        printf("%15s %20.4f \n", "FIX", ffData.angle.Gettheta0(search) *coef);
      else
        printf("%15.6f %20.4f \n", ffData.angle.GetKtheta(search),
               ffData.angle.Gettheta0(search) * coef);

      printed.push_back(angleName);
      printed.push_back(angleNameReverse);
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
    std::string dihedralName, dihedralNameReverse;
    uint dihsize = ffData.dih.GetSizeDih(dName);

    elementNames[0] = kind.atoms[kind.dihedrals[i].a0].type;
    elementNames[1] = kind.atoms[kind.dihedrals[i].a1].type;
    elementNames[2] = kind.atoms[kind.dihedrals[i].a2].type;
    elementNames[3] = kind.atoms[kind.dihedrals[i].a3].type;

    for(uint m = 0; m < ATOMS_PER; ++m) {
      dihedralName.append(elementNames[m]).append("  ");
      dihedralNameReverse.append(elementNames[ATOMS_PER - m - 1]).append("  ");
    }

    if(find(printed.begin(), printed.end(), dihedralName) == printed.end()) {
      for(uint j = 0; j < dihsize; j++) {
        printf("%-20s", dihedralName.c_str());
        printf("%15.6f %4d %15.4f \n", ffData.dih.GetKchi(dName, j),
               ffData.dih.Getn(dName, j),
               ffData.dih.Getdelta(dName, j) * coef);
      }
      printed.push_back(dihedralName);
      printed.push_back(dihedralNameReverse);
      std::cout << endl;
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
      if (i % 20 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.bonds[i].a0 << ' ' << it->second.bonds[i].a1 << ']' << ' ';
    }
    std::cout << std::endl << "\nAngles:";
    for (uint i = 0; i < it->second.angles.size(); i++) {
      if (i % 24 == 0)
        std::cout << std::endl;
      std::cout << "[" << it->second.angles[i].a0 << ' '
                << it->second.angles[i].a1 << ' '
                << it->second.angles[i].a2 << ']' << ' ';
    }
    std::cout << std::endl << "\nDihedrals:";
    for (uint i = 0; i < it->second.dihedrals.size(); i++) {
      if (i % 24 == 0)
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
int ReadPSF(const char* psfFilename, MolMap& kindMap)
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
  ReadPSFAtoms(psf, kindMap, nAtoms);
  //build list of start particles for each type, so we can find it and skip
  //everything else
  std::vector<std::pair<unsigned int, std::string> > firstAtomLookup;
  for (MolMap::iterator it = kindMap.begin();
       it != kindMap.end(); ++it) {
    firstAtomLookup.push_back(std::make_pair(it->second.firstAtomID, it->first));
  }
  std::sort(firstAtomLookup.begin(), firstAtomLookup.end());
  //find bond header+count
  while (strstr(input, "!NBOND") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read bonds from PSF file %s", psfFilename);
      fclose(psf);
      return  READERROR;
    }
  }
  //make sure molecule has bonds, appears before !NBOND
  count = atoi(input);
  if (count != 0) {
    if (ReadPSFBonds(psf, kindMap, firstAtomLookup) == READERROR) {
      fclose(psf);
      return READERROR;
    }
  }
  //find angle header+count
  psf = fopen(psfFilename, "r");
  while (strstr(input, "!NTHETA") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read angles from PSF file %s", psfFilename);
      fclose(psf);
      return READERROR;
    }
  }
  //make sure molecule has angles, count appears before !NTHETA
  count = atoi(input);
  if (count != 0) {
    if (ReadPSFAngles(psf, kindMap, firstAtomLookup) == READERROR) {
      fclose(psf);
      return READERROR;
    }
  }
  //find dihedrals header+count
  psf = fopen(psfFilename, "r");
  while (strstr(input, "!NPHI") == NULL) {
    check = fgets(input, 511, psf);
    if (check == NULL) {
      fprintf(stderr, "ERROR: Unable to read dihedrals from PSF file %s", psfFilename);
      fclose(psf);
      return READERROR;
    }
  }
  //make sure molecule has dihs, count appears before !NPHI
  count = atoi(input);
  if (count != 0) {
    if (ReadPSFDihedrals(psf, kindMap, firstAtomLookup) == READERROR) {
      fclose(psf);
      return READERROR;
    }
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
      it->second.atoms.push_back(Atom(atomName, atomType, charge, mass));
    }
    //still building a molecule...
    else if (it->second.incomplete) {
      if (molID != it->second.firstMolID)
        it->second.incomplete = false;
      else
        it->second.atoms.push_back(Atom(atomName, atomType, charge, mass));
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
                 std::vector<std::pair<unsigned int, std::string> >& firstAtom)
{
  unsigned int atom0, atom1;
  int dummy = fscanf(psf, "%u %u", &atom0, &atom1);
  UNUSED(dummy);
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    //continue if atom has no bonds
    if (currentMol.atoms.size() < 2)
      continue;
    //index of first atom in moleule
    unsigned int molBegin = firstAtom[i].first;
    //index AFTER last atom in molecule
    unsigned int molEnd = molBegin + currentMol.atoms.size();
    while (atom0 < molBegin || atom0 >= molEnd) {
      dummy = fscanf(psf, "%u %u", &atom0, &atom1);
      if (feof(psf) || ferror(psf)) {
        fprintf(stderr, "ERROR: Could not find all bonds in PSF file ");
        return READERROR;
      }
    }
    //read in bonds
    while (atom0 >= molBegin && atom0 < molEnd) {
      currentMol.bonds.push_back(Bond(atom0 - molBegin, atom1 - molBegin));
      dummy = fscanf(psf, "%u %u", &atom0, &atom1);
      if(dummy != 2)
        break;
    }
  }
  return 0;
}

//adds angles in psf to kindMap
//pre: stream is before !NTHETA   post: stream is in angle section just after
//the first appearance of the last molecule
int ReadPSFAngles(FILE* psf, MolMap& kindMap,
                  std::vector<std::pair<unsigned int, std::string> >& firstAtom)
{
  unsigned int atom0, atom1, atom2;
  int dummy = fscanf(psf, "%u %u %u", &atom0, &atom1, &atom2);
  UNUSED(dummy);
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    //continue if atom has no angles
    if (currentMol.atoms.size() < 3)
      continue;
    //index of first atom in moleule
    unsigned int molBegin = firstAtom[i].first;
    //index AFTER last atom in molecule
    unsigned int molEnd = molBegin + currentMol.atoms.size();
    while (atom0 < molBegin || atom0 >= molEnd) {
      dummy = fscanf(psf, "%u %u %u", &atom0, &atom1, &atom2);
      if (feof(psf) || ferror(psf)) {
        fprintf(stderr, "ERROR: Could not find all angles in PSF file ");
        return READERROR;
      }
    }
    //read in angles
    while (atom0 >= molBegin && atom0 < molEnd) {
      currentMol.angles.push_back(Angle(atom0 - molBegin, atom1 - molBegin,
                                        atom2 - molBegin));
      dummy = fscanf(psf, "%u %u %u", &atom0, &atom1, &atom2);
      if(dummy != 3)
        break;
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
                     std::vector<std::pair<unsigned int, std::string> >& firstAtom)
{
  Dihedral dih(0, 0, 0, 0);

  int dummy = fscanf(psf, "%u %u %u %u", &dih.a0, &dih.a1, &dih.a2, &dih.a3);
  UNUSED(dummy);
  //for all atoms
  for (unsigned int i = 0; i < firstAtom.size(); ++i) {
    MolKind& currentMol = kindMap[firstAtom[i].second];
    //continue if molecule has no dihedrals
    if (currentMol.atoms.size() < 4)
      continue;
    //index of first atom in moleule
    unsigned int molBegin = firstAtom[i].first;
    //index AFTER last atom in molecule
    unsigned int molEnd = molBegin + currentMol.atoms.size();
    //continue if molecule has more that 3 atoms but has no dihedrals
    if(i == 0) {
      // if it is the first molecule and index of dihedral is greater than
      // molBegin, it means it does not have any dihedral. It works when
      // we have only two molecule kinds.
      if(dih.a0 > molBegin && dih.a0 > molEnd) {
        continue;
      }
    }
    //scan to to first appearance of molecule
    while (dih.a0 < molBegin || dih.a0 >= molEnd) {
      dummy = fscanf(psf, "%u %u %u %u", &dih.a0, &dih.a1, &dih.a2, &dih.a3);
      if (feof(psf) || ferror(psf)) {
        fprintf(stderr, "ERROR: Could not find all dihedrals in PSF file ");
        return READERROR;
      }
      //for case that molecule has more thatn 3 atoms and represend as
      //second molecule but it has no dihedral. It works when
      // we have only two molecule kinds and no improper
      if (dih.a0 == 0)
        break;
    }
    //read in dihedrals
    while (dih.a0 >= molBegin && dih.a0 < molEnd) {
      dih.a0 -= molBegin;
      dih.a1 -= molBegin;
      dih.a2 -= molBegin;
      dih.a3 -= molBegin;
      //some xplor PSF files have duplicate dihedrals, we need to ignore these
      if (std::find(currentMol.dihedrals.begin(), currentMol.dihedrals.end(),
                    dih) == currentMol.dihedrals.end()) {
        currentMol.dihedrals.push_back(dih);
      }
      dummy = fscanf(psf, "%u %u %u %u", &dih.a0, &dih.a1, &dih.a2, &dih.a3);
      if(dummy != 4)
        break;
    }
  }
  return 0;
}
}
