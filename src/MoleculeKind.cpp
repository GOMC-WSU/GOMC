/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "FFSetup.h"
#include "PDBConst.h" //For resname length.
#include "PRNG.h"
#include "Geometry.h"
#include "Setup.h"
#include "CBMC.h"

#include <vector>
#include <string>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>      //for exit



void MoleculeKind::Init
(std::string const& l_name, Setup const& setup, Forcefield const& forcefield, 
 System& sys)
{
   mol_setup::MolMap::const_iterator dataIterator =
      setup.mol.kindMap.find(l_name);
   if(dataIterator == setup.mol.kindMap.end())
   {
      fprintf(stderr, "Error: Molecule %s not found in PDB file. Exiting.",
              l_name.c_str());
      exit(EXIT_FAILURE);
   }
   const mol_setup::MolKind& molData = dataIterator->second;
   name = l_name;

#if ENSEMBLE == GCMC
   std::map<std::string, double>::const_iterator kindCPIt =
      setup.config.sys.chemPot.cp.find(name);
   chemPot = kindCPIt->second;
#endif
   
   InitAtoms(molData);

   //Once-through topology objects
   nonBonded.Init(molData);
   sortedNB.Init(nonBonded, numAtoms);
   bondList.Init(molData.bonds);
   angles.Init(molData.angles, bondList);
   dihedrals.Init(molData.dihedrals, bondList);

#ifdef VARIABLE_PARTICLE_NUMBER
   builder = cbmc::MakeCBMC(sys, forcefield, *this, setup);
   //builder = new cbmc::LinearVlugt(sys, forcefield, *this, setup);
#endif
}

MoleculeKind::MoleculeKind() : angles(3), dihedrals(4),
   atomMass(NULL), atomCharge(NULL), builder(NULL), 
   atomKind(NULL)
{}

MoleculeKind::~MoleculeKind()
{
   delete[] atomKind;
   delete[] atomMass;
   delete[] atomCharge;
   delete builder;
}




void MoleculeKind::InitAtoms(mol_setup::MolKind const& molData)
{
   numAtoms = molData.atoms.size();
   atomKind = new uint[numAtoms];
   atomMass = new double[numAtoms];
   atomCharge = new double[numAtoms];
   molMass = 0;
   atomNames.clear();

   //convert array of structures to structure of arrays
   for(uint i = 0; i < numAtoms; ++i)
   {
      const mol_setup::Atom& atom = molData.atoms[i];
      atomNames.push_back(atom.name);
      atomMass[i] = atom.mass;
      molMass += atom.mass;
      atomCharge[i] = atom.charge;
      atomKind[i] = atom.kind;
   }
}

