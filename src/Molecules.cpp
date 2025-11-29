/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "Molecules.h"

#include <algorithm> //For count.
#include <cassert>
#include <string>

#include "FFSetup.h" // ""
#include "Forcefield.h"
#include "MolSetup.h" // ""
#include "PDBSetup.h" //For init
#include "Setup.h"
#include "System.h"
#include "VectorLib.h" //for transfer.

class System;

Molecules::Molecules()
    : start(NULL), restartOrderedStart(NULL), kIndex(NULL), countByKind(NULL),
      chain(NULL), kinds(NULL), pairEnCorrections(NULL),
      pairVirCorrections(NULL), printFlag(true) {}

Molecules::~Molecules(void) {
  delete[] start;
  if (restartFromCheckpoint)
    delete[] restartOrderedStart;
  delete[] kIndex;
  delete[] countByKind;
  delete[] chain;
  delete[] beta;
  delete[] occ;
  delete[] kinds;
  delete[] pairEnCorrections;
  delete[] pairVirCorrections;
}

void Molecules::Init(Setup &setup, Forcefield &forcefield, System &sys) {
  pdb_setup::Atoms &atoms = setup.pdb.atoms;
  // Molecule kinds arrays/data.
  kindsCount = setup.mol.kindMap.size();
  countByKind = new uint[kindsCount];
  kinds = new MoleculeKind[kindsCount];
  if (kindsCount != setup.mol.molVars.molKindIndex) {
    std::cout << "Error: Inconsistency between molecule map and number of "
                 "molecule kinds"
              << std::endl
              << "Error: Please report your PDB/PSF files to "
                 "https://github.com/GOMC-WSU/GOMC/issues"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Whether we need to delete the restartOrderedStart array in the destructor
  restartFromCheckpoint = setup.config.in.restart.restartFromCheckpoint;
  // Molecule instance arrays/data
  count = setup.mol.molVars.startIdxMolecules.size();
  atomCount = atoms.beta.size();
  if (count == 0) {
    std::cerr << "Error: No Molecule was found in the PSF file(s)!"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  // chain = new char [atoms.x.size()];
  start = new uint[count + 1];
  start = vect::TransferInto<uint>(start, setup.mol.molVars.startIdxMolecules);
  kIndex = vect::transfer<uint>(setup.mol.molVars.moleculeKinds);
  chain = vect::transfer<char>(atoms.chainLetter);
  beta = vect::transfer<double>(atoms.beta);
  occ = vect::transfer<double>(atoms.occ);

  start[count] = atoms.x.size();
  kIndexCount = setup.mol.molVars.moleculeKinds.size();

  for (uint mk = 0; mk < kindsCount; mk++) {
    countByKind[mk] = std::count(setup.mol.molVars.moleculeNames.begin(),
                                 setup.mol.molVars.moleculeNames.end(),
                                 setup.mol.molVars.moleculeKindNames[mk]);
    kinds[mk].Init(mk, setup.mol.molVars.uniqueMapKeys[mk], setup, forcefield,
                   sys);
  }

#if ENSEMBLE == GCMC
  // check to have all the molecules in psf file that defined in config file
  std::map<std::string, double>::const_iterator
      kindCPIt = setup.config.sys.chemPot.cp.begin(),
      lastOne = setup.config.sys.chemPot.cp.end();

  while (kindCPIt != lastOne) {
    std::string molName = kindCPIt->first;
    mol_setup::MolMap::const_iterator dataIterator = setup.mol.kindMap.begin();
    for (; dataIterator->second.moleculeName != molName &&
           dataIterator != setup.mol.kindMap.end();
         ++dataIterator) {
    }
    if (dataIterator == setup.mol.kindMap.end()) {
      std::cerr << "================================================"
                << std::endl
                << "Error: Molecule " << molName
                << " was not defined in the PDB file." << std::endl
                << "================================================"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    kindCPIt++;
  }
#endif

  if (printFlag) {
    // calculating netcharge of all molecule kind
    double netCharge = 0.0;
    for (uint mk = 0; mk < kindsCount; mk++) {
      netCharge += (countByKind[mk] * kinds[mk].GetMoleculeCharge());
      if (kinds[mk].MoleculeHasCharge()) {
        if (!forcefield.ewald && !forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Ewald Summation method is disabled!\n\n";
        } else if (!forcefield.electrostatic && forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Electrostatic energy calculation is disabled!\n\n";
        }
      }
    }

    if (std::abs(netCharge) > 10E-7) {
      std::cout << "================================================"
                << std::endl
                << std::endl
                << "Warning: Sum of the charge in the system is: " << netCharge
                << std::endl
                << std::endl
                << "================================================"
                << std::endl;
    }
  }

  // Print LJ nonbonded info
  std::vector<uint> totAtomKind;
  std::vector<std::string> atomNames;
  for (uint i = 0; i < kindsCount; ++i) {
    for (uint j = 0; j < kinds[i].NumAtoms(); j++) {
      if (std::find(totAtomKind.begin(), totAtomKind.end(),
                    kinds[i].AtomKind(j)) == totAtomKind.end()) {
        totAtomKind.push_back(kinds[i].AtomKind(j));
        atomNames.push_back(kinds[i].atomTypeNames[j]);
      }
    }
  }

  PrintLJInfo(totAtomKind, atomNames, forcefield);
  printFlag = false;

  // Pair Correction matrices
  pairEnCorrections = new double[kindsCount * kindsCount];
  pairVirCorrections = new double[kindsCount * kindsCount];
  // Initial with zero
  std::fill_n(pairEnCorrections, kindsCount * kindsCount, 0.0);
  std::fill_n(pairVirCorrections, kindsCount * kindsCount, 0.0);

  if (forcefield.useLRC) {
    for (uint i = 0; i < kindsCount; ++i) {
      for (uint j = i; j < kindsCount; ++j) {
        for (uint pI = 0; pI < kinds[i].NumAtoms(); ++pI) {
          for (uint pJ = 0; pJ < kinds[j].NumAtoms(); ++pJ) {
            pairEnCorrections[i * kindsCount + j] +=
                forcefield.particles->EnergyLRC(kinds[i].AtomKind(pI),
                                                kinds[j].AtomKind(pJ));
            pairVirCorrections[i * kindsCount + j] +=
                forcefield.particles->VirialLRC(kinds[i].AtomKind(pI),
                                                kinds[j].AtomKind(pJ));
          }
        }
        // set other side of the diagonal
        pairEnCorrections[j * kindsCount + i] =
            pairEnCorrections[i * kindsCount + j];
        pairVirCorrections[j * kindsCount + i] =
            pairVirCorrections[i * kindsCount + j];
      }
    }
  } else if (forcefield.useIPC) {
    for (uint i = 0; i < kindsCount; ++i) {
      for (uint j = i; j < kindsCount; ++j) {
        for (uint pI = 0; pI < kinds[i].NumAtoms(); ++pI) {
          for (uint pJ = 0; pJ < kinds[j].NumAtoms(); ++pJ) {
            // There is no Impulse energy term. Just Pressure
            pairVirCorrections[i * kindsCount + j] +=
                forcefield.particles->ImpulsePressureCorrection(
                    kinds[i].AtomKind(pI), kinds[j].AtomKind(pJ));
          }
        }
        // set other side of the diagonal
        pairVirCorrections[j * kindsCount + i] =
            pairVirCorrections[i * kindsCount + j];
      }
    }
  }
}

bool Molecules::operator==(const Molecules &other) {
  bool result = true;
  result &= (count == other.count);
  result &= (atomCount == other.atomCount);
  result &= (kindsCount == other.kindsCount);

  for (int m = 0; m < count + 1; ++m) {
    result &= (start[m] == other.start[m]);
  }
  for (int m = 0; m < count; ++m) {
    result &= (kIndex[m] == other.kIndex[m]);
  }
  for (int a = 0; a < atomCount; ++a) {
    result &= (chain[a] == other.chain[a]);
    result &= (beta[a] == other.beta[a]);
    result &= (occ[a] == other.occ[a]);
  }
  for (int k = 0; k < kindsCount; ++k) {
    result &= (countByKind[k] == other.countByKind[k]);
    result &= (kinds[k] == other.kinds[k]);
  }
  for (int k = 0; k < kindsCount * kindsCount; ++k) {
    result &= (pairEnCorrections[k] == other.pairEnCorrections[k]);
    result &= (pairVirCorrections[k] == other.pairVirCorrections[k]);
  }
  return result;
}

void Molecules::PrintLJInfo(std::vector<uint> &totAtomKind,
                            std::vector<std::string> &names,
                            Forcefield &forcefield) {
  if (printFlag) {
    uint size = totAtomKind.size();
    printf("NonBonded 1-4 parameters:\n");
    if (forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for (uint i = 0; i < size; i++) {
      for (uint j = i; j < size; j++) {
        if (forcefield.exp6) {
          printf(
              "%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
              names[i].c_str(), names[j].c_str(),
              forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                                                   totAtomKind[j]),
              forcefield.particles->GetSigma_1_4(totAtomKind[i],
                                                 totAtomKind[j]),
              forcefield.particles->GetRmax_1_4(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetRmin_1_4(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetN_1_4(totAtomKind[i], totAtomKind[j]));
        } else {
          printf(
              "%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
              names[j].c_str(),
              forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                                                   totAtomKind[j]),
              forcefield.particles->GetSigma_1_4(totAtomKind[i],
                                                 totAtomKind[j]),
              forcefield.particles->GetN_1_4(totAtomKind[i], totAtomKind[j]));
        }
      }
    }

    std::cout << std::endl;
    printf("NonBonded parameters:\n");
    if (forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for (uint i = 0; i < size; i++) {
      for (uint j = i; j < size; j++) {
        if (forcefield.exp6) {
          printf(
              "%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
              names[i].c_str(), names[j].c_str(),
              forcefield.particles->GetEpsilon(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetSigma(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetRmax(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetRmin(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetN(totAtomKind[i], totAtomKind[j]));
        } else {
          printf(
              "%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
              names[j].c_str(),
              forcefield.particles->GetEpsilon(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetSigma(totAtomKind[i], totAtomKind[j]),
              forcefield.particles->GetN(totAtomKind[i], totAtomKind[j]));
        }
      }
    }
    std::cout << std::endl;
  }
}
