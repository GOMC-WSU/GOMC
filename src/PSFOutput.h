/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef PSF_OUTPUT_H
#define PSF_OUTPUT_H

#include <string>
#include <vector>

#include "AlphaNum.h"
#include "BasicTypes.h"
#include "MolSetup.h"
#include "MoleculeLookup.h"
#include "OutputAbstracts.h"
#include "System.h"
#include "Writer.h"

class Molecules;

class PSFOutput : public OutputableBase {
public:
  PSFOutput(const Molecules &molecules, const System &sys, Setup &set);

  // Output PSF file to filename using default remarks
  void PrintPSF(const std::string &filename) const;

  // Output PSF file to filename, specifying remarks
  void PrintPSF(const std::string &filename,
                const std::vector<std::string> &remarks) const;

  // PSF does not need to sample on every step, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  const Molecules *molecules;
  const MoleculeLookup &molLookRef;
  std::vector<mol_setup::MolKind> molKinds;
  std::vector<std::string> molNames;
  uint totalAtoms;
  uint totalBonds;
  uint totalAngles;
  uint totalDihs;
  uint totalImps;
  uint totalDons;
  uint totalAccs;
  uint totalNNBs;
  uint totalGrps;
  uint totalCrtrms;

  uint boxAtoms[BOX_TOTAL];
  uint boxBonds[BOX_TOTAL];
  uint boxAngles[BOX_TOTAL];
  uint boxDihs[BOX_TOTAL];
  uint boxImps[BOX_TOTAL];
  uint boxDons[BOX_TOTAL];
  uint boxAccs[BOX_TOTAL];
  uint boxNNBs[BOX_TOTAL];
  uint boxGrps[BOX_TOTAL];
  uint boxCrtrms[BOX_TOTAL];

  void PrintRemarks(FILE *outfile,
                    const std::vector<std::string> &remarks) const;
  void PrintAtoms(FILE *outfile) const;
  void PrintBonds(FILE *outfile) const;
  void PrintAngles(FILE *outfile) const;
  void PrintDihedrals(FILE *outfile) const;
  void PrintImpropers(FILE *outfile) const;
  void PrintDonors(FILE *outfile) const;
  void PrintAcceptors(FILE *outfile) const;
  void PrintExplicitNonbondedExclusions(FILE *outfile) const;
  void PrintGroups(FILE *outfile) const;
  void PrintCrossTerms(FILE *outfile) const;

  void PrintRemarksInBox(FILE *outfile, uint b) const;
  void PrintAtomsInBox(FILE *outfile, uint b) const;
  void PrintBondsInBox(FILE *outfile, uint b) const;
  void PrintAnglesInBox(FILE *outfile, uint b) const;
  void PrintDihedralsInBox(FILE *outfile, uint b) const;
  void PrintImpropersInBox(FILE *outfile, uint b) const;
  void PrintDonorsInBox(FILE *outfile, uint b) const;
  void PrintAcceptorsInBox(FILE *outfile, uint b) const;
  void PrintExplicitNonbondedExclusionsInBox(FILE *outfile, uint b) const;
  void PrintGroupsInBox(FILE *outfile, uint b) const;
  void PrintCrossTermsInBox(FILE *outfile, uint b) const;

  void PrintNAMDCompliantSuffix(FILE *outfile) const;
  void PrintNAMDCompliantSuffixInBox(FILE *outfile) const;

  void CountMolecules();
  void CountMoleculesInBoxes();

  FILE *outRebuildRestart[BOX_TOTAL];
  std::string outFName;
  std::string outRebuildRestartFName[BOX_TOTAL];
  std::vector<std::string> moleculeSegmentNames;
};

#endif /*PSF_OUTPUT_H*/
