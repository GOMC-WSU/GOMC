/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Molecules.h"
#include "Setup.h"
#include "PDBSetup.h" //For init
#include "MolSetup.h" // ""
#include "FFSetup.h"  // ""
#include "Forcefield.h"
#include "VectorLib.h" //for transfer.
#include <algorithm> //For count.
#include <string>
#include "System.h"

class System;


Molecules::Molecules() : start(NULL), kIndex(NULL), countByKind(NULL),
  chain(NULL), kinds(NULL), pairEnCorrections(NULL),
  pairVirCorrections(NULL), printFlag(true), sortedMoleculeSegmentName(NULL), 
  sortedMoleculeIndices(NULL){}

Molecules::~Molecules(void)
{
  delete[] start;
  delete[] kIndex;
  delete[] countByKind;
  delete[] chain;
  delete[] kinds;
  delete[] pairEnCorrections;
  delete[] pairVirCorrections;
  delete[] sortedMoleculeSegmentName;
  delete[] sortedMoleculeIndices;
}

void Molecules::Init(Setup & setup, Forcefield & forcefield,
                     System & sys)
{
  pdb_setup::Atoms& atoms = setup.pdb.atoms;
  //Molecule kinds arrays/data.
  kindsCount = setup.mol.kindMap.size();
  countByKind = new uint[kindsCount];
  kinds = new MoleculeKind[kindsCount];

  //Molecule instance arrays/data
  count = setup.mol.molVars.startIdxMolecules.size();
  if (count == 0) {
    std::cerr << "Error: No Molecule was found in the PSF file(s)!" << std::endl;
    exit(EXIT_FAILURE);
  }
  //chain = new char [atoms.x.size()];
  start = new uint [count + 1];
  sortedMoleculeIndices = new uint [count];

  enableGenerateSegmentOut = setup.mol.molVars.enableGenerateSegmentOut;
  enableSortedSegmentOut = setup.mol.molVars.enableSortedSegmentOut;

  if(enableSortedSegmentOut){
    /* We need to create the sortedArray in the method and return it, instead of directly modifying the class' vector
     because of the way this class is initialized twice */
    SortMoleculesBySegment(setup.mol.molVars.moleculeSegmentNames);
  }

  start = vect::TransferInto<uint>(start, setup.mol.molVars.startIdxMolecules);
  kIndex = vect::transfer<uint>(setup.mol.molVars.moleculeKinds);
  chain = vect::transfer<char>(atoms.chainLetter);

  start[count] = atoms.x.size();
  kIndexCount = setup.mol.molVars.moleculeKinds.size();

  for (uint mk = 0 ; mk < kindsCount; mk++) {
    countByKind[mk] =
      std::count(setup.mol.molVars.moleculeNames.begin(), setup.mol.molVars.moleculeNames.end(),
                 setup.mol.molVars.moleculeKindNames[mk]);
    kinds[mk].Init(setup.mol.molVars.moleculeKindNames[mk], setup, forcefield, sys);
  }

#if ENSEMBLE == GCMC
  //check to have all the molecules in psf file that defined in config file
  std::map<std::string, double>::const_iterator kindCPIt =
    setup.config.sys.chemPot.cp.begin(),
    lastOne = setup.config.sys.chemPot.cp.end();

  while(kindCPIt != lastOne) {
    std::string molName = kindCPIt->first;
    mol_setup::MolMap::const_iterator dataIterator =
      setup.mol.kindMap.find(molName);
    if(dataIterator == setup.mol.kindMap.end()) {
      std::cerr << "================================================"
                << std::endl << "Error: Molecule " << molName
                << " was not defined in the PDB file." << std::endl
                << "================================================"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    kindCPIt++;
  }
#endif

  if(printFlag) {
    //calculating netcharge of all molecule kind
    double netCharge = 0.0;
    for (uint mk = 0 ; mk < kindsCount; mk++) {
      netCharge += (countByKind[mk] * kinds[mk].GetMoleculeCharge());
      if(kinds[mk].MoleculeHasCharge()) {
        if(!forcefield.ewald && !forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Ewald Summation method is disabled!\n\n";
        } else if(!forcefield.electrostatic && forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Electrostatic energy calculation is disabled!\n\n";
        }
      }
    }

    if(std::abs(netCharge) > 10E-7) {
      std::cout << "================================================"
                << std::endl << std::endl
                << "Warning: Sum of the charge in the system is: "
                << netCharge << std::endl << std::endl
                <<  "================================================"
                << std::endl;

    }
  }

  //Print LJ nonbonded info
  std::vector<uint> totAtomKind;
  std::vector<std::string> atomNames;
  for(uint i = 0; i < kindsCount; ++i) {
    for(uint j = 0; j < kinds[i].NumAtoms(); j++) {
      if(std::find(totAtomKind.begin(), totAtomKind.end(), kinds[i].AtomKind(j))
          == totAtomKind.end()) {
        totAtomKind.push_back(kinds[i].AtomKind(j));
        atomNames.push_back(kinds[i].atomTypeNames[j]);
      }
    }
  }

  PrintLJInfo(totAtomKind, atomNames, forcefield);
  printFlag = false;

  //Pair Correction matrices
  pairEnCorrections = new double[kindsCount * kindsCount];
  pairVirCorrections = new double[kindsCount * kindsCount];
  for(uint i = 0; i < kindsCount; ++i) {
    for(uint j = i; j < kindsCount; ++j) {
      pairEnCorrections[i * kindsCount + j] = 0.0;
      pairVirCorrections[i * kindsCount + j] = 0.0;
      for(uint pI = 0; pI < kinds[i].NumAtoms(); ++pI) {
        for(uint pJ = 0; pJ < kinds[j].NumAtoms(); ++pJ) {
          pairEnCorrections[i * kindsCount + j] +=
            forcefield.particles->EnergyLRC(kinds[i].AtomKind(pI),
                                            kinds[j].AtomKind(pJ));
          pairVirCorrections[i * kindsCount + j] +=
            forcefield.particles->VirialLRC(kinds[i].AtomKind(pI),
                                            kinds[j].AtomKind(pJ));
        }
      }
      //set other side of the diagonal
      pairEnCorrections[j * kindsCount + i] =
        pairEnCorrections[i * kindsCount + j];
      pairVirCorrections[j * kindsCount + i] =
        pairVirCorrections[i * kindsCount + j];
    }
  }
}

void Molecules::PrintLJInfo(std::vector<uint> &totAtomKind,
                            std::vector<std::string> &names,
                            Forcefield & forcefield)
{
  if(printFlag) {
    uint size =  totAtomKind.size();
    printf("NonBonded 1-4 parameters:\n");
    if(forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for(uint i = 0; i < size; i++) {
      for(uint j = i; j < size; j++) {
        if(forcefield.exp6) {
          printf("%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
                 names[i].c_str(), names[j].c_str(),
                 forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetSigma_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetRmax_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetRmin_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetN_1_4(totAtomKind[i],
                                                totAtomKind[j]));
        } else {
          printf("%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
                 names[j].c_str(),
                 forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetSigma_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetN_1_4(totAtomKind[i],
                                                totAtomKind[j]));
        }
      }
    }

    std::cout << std::endl;
    printf("NonBonded parameters:\n");
    if(forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for(uint i = 0; i < size; i++) {
      for(uint j = i; j < size; j++) {
        if(forcefield.exp6) {
          printf("%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
                 names[i].c_str(), names[j].c_str(),
                 forcefield.particles->GetEpsilon(totAtomKind[i],
                                                  totAtomKind[j]),
                 forcefield.particles->GetSigma(totAtomKind[i],
                                                totAtomKind[j]),
                 forcefield.particles->GetRmax(totAtomKind[i],
                                               totAtomKind[j]),
                 forcefield.particles->GetRmin(totAtomKind[i],
                                               totAtomKind[j]),
                 forcefield.particles->GetN(totAtomKind[i],
                                            totAtomKind[j]));
        } else {
          printf("%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
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

void  Molecules::SortMoleculesBySegment(std::vector<std::string> & unorderedSegments){

  /* For Hybid MC-MD Cycle Consistency between molecular order
     Sort these three vectors according to alphanmeric segment label */ 

  std::vector<std::string> sortedSegmentName;
  std::vector<uint> sortedSegmentIndices;

  std::vector<uint> unsortedSegmentIndices(count);
  std::iota(unsortedSegmentIndices.begin(), unsortedSegmentIndices.end(), 0);
  //declaring vector of pairs
  std::vector< std::pair <std::string,int> > pairVector;
  // Entering values in vector of pairs
  for (int i=0; i<count; i++)
      pairVector.push_back( std::make_pair(unorderedSegments[i],unsortedSegmentIndices[i]) );

  // Using simple sort() function to sort
  std::sort(pairVector.begin(), pairVector.end(), compare());

  for (int i = 0; i < count; i++){
    sortedSegmentName.push_back(unorderedSegments[pairVector[i].second]); 
    sortedSegmentIndices.push_back(pairVector[i].second); 
  }


  vect::TransferInto<uint>(this->sortedMoleculeIndices, sortedSegmentIndices);
  //chain = vect::transfer<char>(sortedChain);
  //return sortedSegmentIndices;
  /* For Hybid MC-MD Cycle Consistency between molecular order
    Sort these three vectors according to alphanmeric segment label */ 

}

