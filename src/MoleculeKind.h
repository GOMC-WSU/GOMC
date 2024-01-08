/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOLECULE_KIND_H
#define MOLECULE_KIND_H

#include <cassert>
#include <string>
#include <vector>

#include "BasicTypes.h" //For uint
#include "CBMC.h"
#include "EnsemblePreprocessor.h"
#include "Geometry.h" //members
#include "SubdividedArray.h"

namespace mol_setup {
class Atom;
class MolKind;
} // namespace mol_setup
namespace ff_setup {
class Bond;
class FFBase;
} // namespace ff_setup
namespace cbmc {
class TrialMol;
}

class FFSetup;
class PRNG;
struct MolPick;
class System;
class Forcefield;
class Setup;

class MoleculeKind {
public:
  MoleculeKind();
  ~MoleculeKind();

  uint NumAtoms() const { return numAtoms; }
  uint NumBonds() const { return bondList.count; }
  uint NumAngles() const { return angles.Count(); }
  uint NumDihs() const { return dihedrals.Count(); }
  uint NumImps() const { return impropers.Count(); }
  uint NumDons() const { return donorList.count; }
  uint NumAccs() const { return acceptorList.count; }
  /*
  uint NumNNBs() const
  {
    return nnbs.Count();
  }
  uint NumGrps() const
  {
    return groups.Count();
  }
  uint NumCrtrms() const
  {
    return crossterms.Count();
  }
  */
  uint AtomKind(const uint a) const { return atomKind[a]; }
  double AtomCharge(const uint a) const { return atomCharge[a]; }

  // Initialize this kind
  // Exits program if param and psf files don't match
  void Init(uint &l_kindIndex, std::string const &l_name, Setup const &setup,
            Forcefield const &forcefield, System &sys);

  // Invoke CBMC, oldMol and newMol will be modified
  void Build(cbmc::TrialMol &oldMol, cbmc::TrialMol &newMol,
             const uint molIndex) {
    builder->Build(oldMol, newMol, molIndex);
  }

  // CBMC for regrowth move
  void Regrowth(cbmc::TrialMol &oldMol, cbmc::TrialMol &newMol,
                const uint molIndex) {
    builder->Regrowth(oldMol, newMol, molIndex);
  }

  // Crank Shaft move
  void CrankShaft(cbmc::TrialMol &oldMol, cbmc::TrialMol &newMol,
                  const uint molIndex) {
    builder->CrankShaft(oldMol, newMol, molIndex);
  }

  // Targeted Swap move
  void BuildGrowInCav(cbmc::TrialMol &oldMol, cbmc::TrialMol &newMol,
                      const uint molIndex) {
    builder->BuildGrowInCav(oldMol, newMol, molIndex);
  }

  // Used in MEMC move
  void BuildIDNew(cbmc::TrialMol &newMol, const uint molIndex) {
    builder->BuildIDNew(newMol, molIndex);
  }

  void BuildIDOld(cbmc::TrialMol &oldMol, const uint molIndex) {
    builder->BuildIDOld(oldMol, molIndex);
  }

  void BuildNew(cbmc::TrialMol &newMol, const uint molIndex) {
    builder->BuildNew(newMol, molIndex);
  }

  void BuildOld(cbmc::TrialMol &oldMol, const uint molIndex) {
    builder->BuildOld(oldMol, molIndex);
  }

  void BuildGrowNew(cbmc::TrialMol &newMol, const uint molIndex) {
    builder->BuildGrowNew(newMol, molIndex);
  }

  void BuildGrowOld(cbmc::TrialMol &oldMol, const uint molIndex) {
    builder->BuildGrowOld(oldMol, molIndex);
  }

  double GetMoleculeCharge();

  bool MoleculeHasCharge();

  bool operator==(const MoleculeKind &other);

  SortedNonbond sortedNB, sortedNB_1_4, sortedNB_1_3, sortedEwaldNB;

  // these are used for total energy calculations, see Geometry.h/cpp
  Nonbond nonBonded;
  Nonbond_1_4 nonBonded_1_4;
  Nonbond_1_3 nonBonded_1_3;
  EwaldNonbond nonEwaldBonded;

  BondList bondList, donorList, acceptorList;
  GeomFeature angles;
  GeomFeature dihedrals;
  GeomFeature impropers;

  bool oneThree, oneFour;

  // uniqueName - guarunteed to be unique, the map key
  // name - not guarunteed to be unique
  //  -if a protein, == PROT(A...Z)
  //  -if non-protein, residue name
  std::string name, uniqueName;
  ;
  uint kindIndex;
  std::vector<std::string> atomNames, atomTypeNames, resNames;
  double molMass;

  double *atomMass;

  bool isMultiResidue;
  std::vector<uint> intraMoleculeResIDs;

#if ENSEMBLE == GCMC
  double chemPot;
#endif

private:
  void InitAtoms(mol_setup::MolKind const &molData);

  // uses buildBonds to check if molecule is branched
  // bool CheckBranches();
  void InitCBMC(System &sys, Forcefield &ff, Setup &set);

  cbmc::CBMC *builder;

  uint numAtoms;
  uint *atomKind;
  double *atomCharge;
};

#endif /*MOLECULE_KIND_H*/
