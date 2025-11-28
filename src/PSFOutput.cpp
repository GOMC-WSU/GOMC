/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "PSFOutput.h"

#include <cstdio>

#include "Molecules.h"

using namespace mol_setup;

namespace {
const char *remarkHeader = "!NTITLE";
const char *atomHeader = "!NATOM";
const char *bondHeader = "!NBOND: bonds";
const char *angleHeader = "!NTHETA: angles";
const char *dihedralHeader = "!NPHI: dihedrals";
const char *improperHeader = "!NIMPHI: impropers";
const char *donorHeader = "!NDON: donors";
const char *acceptorHeader = "!NACC: acceptors";
const char *excludedHeader = "!NNB";
const char *groupHeader = "!NGRP";

const char *headerFormat = "%8d %s \n";
// atom ID, segment name, residue ID, residue name,
// atom name, atom type, charge, mass, and an unused 0
// const char* atomFormat = "%8d%4s%3d%7s%4s%6s%12.6f%14.4f%12d\n";
const char *atomFormat = "%8d %-5s%-5d%-5s%-5s%-3s%12.6f%14.4f%12d\n";

const int bondPerLine = 4;
const int anglePerLine = 3;
const int dihPerLine = 2;
const int impsPerLine = 2;
const int donorPerLine = 4;
const int acceptorPerLine = 4;

} // namespace

PSFOutput::PSFOutput(const Molecules &molecules, const System &sys, Setup &set)
    : molecules(&molecules), molLookRef(sys.molLookup),
      molNames(set.mol.molVars.moleculeKindNames),
      moleculeSegmentNames(set.mol.molVars.moleculeSegmentNames) {
  molKinds.resize(set.mol.kindMap.size());
  for (uint i = 0; i < set.mol.molVars.uniqueMapKeys.size(); ++i) {
    molKinds[i] = set.mol.kindMap[set.mol.molVars.uniqueMapKeys[i]];
  }
  outFName = set.config.out.state.files.psf.name;
  /* To eliminate arithmetic exceptions */
  stepsPerOut = 1;
  printOnFirstStep = true;
}

void PSFOutput::Init(pdb_setup::Atoms const &atoms,
                     config_setup::Output const &output) {
  std::string bStr = "", aliasStr = "", numStr = "";
  sstrm::Converter toStr;
  enableRestOut = output.restart.settings.enable || forceOutput;
  stepsRestPerOut = output.restart.settings.frequency;
  if (enableRestOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      // Get alias string, based on box #.
      bStr = "Box ";
      numStr = "";
      toStr << b + 1;
      toStr >> numStr;
      aliasStr = "Output PSF file for Box ";
      aliasStr += numStr;
      // NEW_RESTART_COD
      outRebuildRestartFName[b] = output.state.files.splitPSF.name[b];
      std::string newStrAddOn = "_restart.psf";
      outRebuildRestartFName[b].replace(outRebuildRestartFName[b].end() - 4,
                                        outRebuildRestartFName[b].end(),
                                        newStrAddOn);
    }
  }
}

void PSFOutput::DoOutputRestart(const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::PSF_RESTART_OUTPUT);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    FILE *outfile = fopen(outRebuildRestartFName[b].c_str(), "w");
    if (outfile == NULL) {
      fprintf(stderr, "Error opening PSF output file %s",
              outRebuildRestartFName[b].c_str());
      return;
    }
    fprintf(outfile, "PSF\n\n");
    CountMoleculesInBoxes();
    PrintRemarksInBox(outfile, b);
    PrintAtomsInBox(outfile, b);
    PrintBondsInBox(outfile, b);
    PrintAnglesInBox(outfile, b);
    PrintDihedralsInBox(outfile, b);
    PrintImpropersInBox(outfile, b);
    PrintDonorsInBox(outfile, b);
    PrintAcceptorsInBox(outfile, b);
    PrintNAMDCompliantSuffixInBox(outfile);
    fclose(outfile);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::PSF_RESTART_OUTPUT);
}

/* Output (merged_psf) occurs in Constructor only */
void PSFOutput::DoOutput(const ulong step) {
  /* merged psf only prints on first step
     Don't print on Checkpoint restarts */
  if (restartFromCheckpoint || step != startStep)
    return;

  GOMC_EVENT_START(1, GomcProfileEvent::PSF_MERGED_OUTPUT);

  CountMolecules();

  std::vector<std::string> remarks;
  // default file remarks
  remarks.push_back("Combined PSF produced by GOMC");
  remarks.push_back(
      "Contains Geometry data for molecules in ALL boxes in the system");
  FILE *outfile = fopen(outFName.c_str(), "w");
  if (outfile == NULL) {
    fprintf(stderr, "Error opening PSF output file %s", outFName.c_str());
    return;
  }

  fprintf(outfile, "PSF\n\n");
  PrintRemarks(outfile, remarks);
  PrintAtoms(outfile);
  PrintBonds(outfile);
  PrintAngles(outfile);
  PrintDihedrals(outfile);
  PrintImpropers(outfile);
  PrintDonors(outfile);
  PrintAcceptors(outfile);
  PrintNAMDCompliantSuffix(outfile);
  fclose(outfile);

  std::cout << "Printed combined psf to file " << outFName << '\n';

  GOMC_EVENT_STOP(1, GomcProfileEvent::PSF_MERGED_OUTPUT);
}

void PSFOutput::CountMolecules() {
  totalAngles = 0;
  totalAtoms = 0;
  totalBonds = 0;
  totalDihs = 0;
  totalImps = 0;
  totalDons = 0;
  totalAccs = 0;
  uint atomT = 0;

  for (uint b = 0; b < BOX_TOTAL; b++) {
    for (uint k = 0; k < molKinds.size(); ++k) {
      // This doesnt work when the molecules are interspersed instead of all of
      // one type then all the other
      // const MoleculeKind& molKind = molecules->GetKind(atomT);

      const MoleculeKind &molKind = molecules->kinds[k];

      totalAtoms += molKind.NumAtoms() * molLookRef.NumKindInBox(k, b);
      totalBonds += molKind.NumBonds() * molLookRef.NumKindInBox(k, b);
      totalAngles += molKind.NumAngles() * molLookRef.NumKindInBox(k, b);
      totalDihs += molKind.NumDihs() * molLookRef.NumKindInBox(k, b);
      totalImps += molKind.NumImps() * molLookRef.NumKindInBox(k, b);
      totalDons += molKind.NumDons() * molLookRef.NumKindInBox(k, b);
      totalAccs += molKind.NumAccs() * molLookRef.NumKindInBox(k, b);
      /*
      totalNNBs += molKind.NumNNBs() * molLookRef.NumKindInBox(k, b);
      totalGrps += molKind.NumGrps() * molLookRef.NumKindInBox(k, b);
      totalCrtrms += molKind.NumCrtrms() * molLookRef.NumKindInBox(k, b);
      */
      atomT += molLookRef.NumKindInBox(k, b);
    }
  }
}

void PSFOutput::CountMoleculesInBoxes() {
  uint atomT = 0;

  for (uint b = 0; b < BOX_TOTAL; b++) {
    boxAtoms[b] = 0;
    boxBonds[b] = 0;
    boxAngles[b] = 0;
    boxDihs[b] = 0;
    boxImps[b] = 0;
    boxDons[b] = 0;
    boxAccs[b] = 0;
    boxNNBs[b] = 0;
    boxGrps[b] = 0;
    boxCrtrms[b] = 0;
    for (uint k = 0; k < molKinds.size(); ++k) {
      // This doesnt work when the molecules are interspersed instead of all of
      // one type then all the other
      // const MoleculeKind& molKind = molecules->GetKind(atomT);
      const MoleculeKind &molKind = molecules->kinds[k];

      boxAtoms[b] += molKind.NumAtoms() * molLookRef.NumKindInBox(k, b);
      boxBonds[b] += molKind.NumBonds() * molLookRef.NumKindInBox(k, b);
      boxAngles[b] += molKind.NumAngles() * molLookRef.NumKindInBox(k, b);
      boxDihs[b] += molKind.NumDihs() * molLookRef.NumKindInBox(k, b);
      boxImps[b] += molKind.NumImps() * molLookRef.NumKindInBox(k, b);
      boxDons[b] += molKind.NumDons() * molLookRef.NumKindInBox(k, b);
      boxAccs[b] += molKind.NumAccs() * molLookRef.NumKindInBox(k, b);
      /*
      boxNNBs[b] += molKind.NumNNBs() * molLookRef.NumKindInBox(k, b);
      boxGrps[b] += molKind.NumGrps() * molLookRef.NumKindInBox(k, b);
      boxCrtrms[b] += molKind.NumCrtrms() * molLookRef.NumKindInBox(k, b);
      */
    }
  }
}

void PSFOutput::PrintPSF(const std::string &filename) const {
  GOMC_EVENT_START(1, GomcProfileEvent::PSF_MERGED_OUTPUT);
  std::vector<std::string> remarks;
  // default file remarks
  remarks.push_back("Combined PSF produced by GOMC");
  remarks.push_back(
      "Contains Geometry data for molecules in ALL boxes in the system");
  PrintPSF(filename, remarks);
  GOMC_EVENT_STOP(1, GomcProfileEvent::PSF_MERGED_OUTPUT);
}

void PSFOutput::PrintPSF(const std::string &filename,
                         const std::vector<std::string> &remarks) const {
  FILE *outfile = fopen(filename.c_str(), "w");
  if (outfile == NULL) {
    fprintf(stderr, "Error opening PSF output file %s", filename.c_str());
    return;
  }

  fprintf(outfile, "PSF\n\n");
  PrintRemarks(outfile, remarks);
  PrintAtoms(outfile);
  PrintBonds(outfile);
  PrintAngles(outfile);
  PrintDihedrals(outfile);
  /* Imps, Dons, Accs are simply read and printed,
    They don't influence a GOMC simulation */
  PrintImpropers(outfile);
  PrintDonors(outfile);
  PrintAcceptors(outfile);
  PrintNAMDCompliantSuffix(outfile);
  fclose(outfile);
}

void PSFOutput::PrintRemarks(FILE *outfile,
                             const std::vector<std::string> &remarks) const {
  fprintf(outfile, headerFormat, remarks.size(), remarkHeader);
  for (uint i = 0; i < remarks.size(); ++i) {
    fprintf(outfile, " REMARKS ");
    fprintf(outfile, "%s", remarks[i].c_str());
    fputc('\n', outfile);
  }
  fputc('\n', outfile);
}

void PSFOutput::PrintAtoms(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalAtoms, atomHeader);
  // silly psfs index from 1
  uint atomID = 1;
  uint resID = 1;
  uint thisKIndex = 0, nAtoms = 0, mI = 0;
  uint pStart = 0, pEnd = 0;
  // Start particle numbering @ 1
  for (uint mol = 0; mol < molecules->count; ++mol) {
    // If this isn't checkpoint restarted, then this is
    thisKIndex = molecules->kIndex[mol];
    nAtoms = molKinds[thisKIndex].atoms.size();

    for (uint at = 0; at < nAtoms; ++at) {
      const Atom *thisAtom = &molKinds[thisKIndex].atoms[at];
      // atom ID, segment name, residue ID, residue name,
      // atom name, atom type, charge, mass, and an unused 0

      if (molKinds[thisKIndex].isMultiResidue) {
        fprintf(outfile, atomFormat, atomID, moleculeSegmentNames[mol].c_str(),
                resID + molKinds[thisKIndex].intraMoleculeResIDs[at],
                thisAtom->residue.c_str(), thisAtom->name.c_str(),
                thisAtom->type.c_str(), thisAtom->charge, thisAtom->mass, 0);
      } else {
        fprintf(outfile, atomFormat, atomID, moleculeSegmentNames[mol].c_str(),
                resID, thisAtom->residue.c_str(), thisAtom->name.c_str(),
                thisAtom->type.c_str(), thisAtom->charge, thisAtom->mass, 0);
      }

      ++atomID;
    }
    /* This isn't actually residue, it is running count of the number of
      molecule kinds we have printed */
    ++resID;
    /* To add additional intramolecular residues */
    if (molKinds[thisKIndex].isMultiResidue) {
      resID += molKinds[thisKIndex].intraMoleculeResIDs.back();
    }

    if (resID == 10000)
      resID = 1;
  }
  fputc('\n', outfile);
}

void PSFOutput::PrintBonds(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalBonds, bondHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    // If this isn't checkpoint restarted, then this is
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.bonds.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.bonds[i].a0 + atomID,
              thisKind.bonds[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == bondPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintAngles(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalAngles, angleHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    // If this isn't checkpoint restarted, then this is
    // mI = *m;
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.angles.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d", thisKind.angles[i].a0 + atomID,
              thisKind.angles[i].a1 + atomID, thisKind.angles[i].a2 + atomID);
      ++lineEntry;
      if (lineEntry == anglePerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}
void PSFOutput::PrintDihedrals(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalDihs, dihedralHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.dihedrals.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d%8d", thisKind.dihedrals[i].a0 + atomID,
              thisKind.dihedrals[i].a1 + atomID,
              thisKind.dihedrals[i].a2 + atomID,
              thisKind.dihedrals[i].a3 + atomID);
      ++lineEntry;
      if (lineEntry == dihPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintImpropers(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalImps, improperHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.impropers.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d%8d", thisKind.impropers[i].a0 + atomID,
              thisKind.impropers[i].a1 + atomID,
              thisKind.impropers[i].a2 + atomID,
              thisKind.impropers[i].a3 + atomID);
      ++lineEntry;
      if (lineEntry == impsPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintDonors(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalDons, donorHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    // If this isn't checkpoint restarted, then this is
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.donors.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.donors[i].a0 + atomID,
              thisKind.donors[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == donorPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintAcceptors(FILE *outfile) const {
  fprintf(outfile, headerFormat, totalAccs, acceptorHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  uint thisKIndex = 0, mI = 0;
  for (uint mol = 0; mol < molecules->count; ++mol) {
    // If this isn't checkpoint restarted, then this is
    thisKIndex = molecules->kIndex[mol];
    const MolKind &thisKind = molKinds[thisKIndex];
    for (uint i = 0; i < thisKind.acceptors.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.acceptors[i].a0 + atomID,
              thisKind.acceptors[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == acceptorPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintNAMDCompliantSuffix(FILE *outfile) const {
  // fprintf(outfile, headerFormat, 0, improperHeader);
  // fputs("\n\n", outfile);
  // fprintf(outfile, headerFormat, 0, donorHeader);
  // fputs("\n\n", outfile);
  // fprintf(outfile, headerFormat, 0, acceptorHeader);
  // fputs("\n\n", outfile);
  fprintf(outfile, headerFormat, 0, excludedHeader);
  fputs("\n\n", outfile);
  fprintf(outfile, headerFormat, 0, groupHeader);
}

void PSFOutput::PrintRemarksInBox(FILE *outfile, uint b) const {
  std::vector<std::string> remarks;
  std::string boxSpecific;
  // default file remarks
  remarks.push_back("Restart PSF produced by GOMC");
  boxSpecific = std::string("Contains Geometry data for molecules in Box " +
                            std::to_string(b));
  remarks.push_back(boxSpecific);
  PrintRemarks(outfile, remarks);
}
void PSFOutput::PrintAtomsInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxAtoms[b], atomHeader);
  // silly psfs index from 1
  uint atomID = 1;
  uint resID = 1;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    uint thisKind = molecules->kIndex[*thisMol];
    uint nAtoms = molKinds[thisKind].atoms.size();

    for (uint at = 0; at < nAtoms; ++at) {
      const Atom *thisAtom = &molKinds[thisKind].atoms[at];
      // atom ID, segment name, residue ID, residue name,
      // atom name, atom type, charge, mass, and an unused 0

      if (molKinds[thisKind].isMultiResidue) {
        fprintf(outfile, atomFormat, atomID,
                moleculeSegmentNames[*thisMol].c_str(),
                resID + molKinds[thisKind].intraMoleculeResIDs[at],
                thisAtom->residue.c_str(), thisAtom->name.c_str(),
                thisAtom->type.c_str(), thisAtom->charge, thisAtom->mass, 0);
      } else {
        fprintf(outfile, atomFormat, atomID,
                moleculeSegmentNames[*thisMol].c_str(), resID,
                thisAtom->residue.c_str(), thisAtom->name.c_str(),
                thisAtom->type.c_str(), thisAtom->charge, thisAtom->mass, 0);
      }
      ++atomID;
    }
    ++resID;
    /* To add additional intramolecular residues */
    if (molKinds[thisKind].isMultiResidue) {
      resID += molKinds[thisKind].intraMoleculeResIDs.back();
    }

    // ???
    if (resID == 10000)
      resID = 1;
  }
  fputc('\n', outfile);
}
void PSFOutput::PrintBondsInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxBonds[b], bondHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.bonds.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.bonds[i].a0 + atomID,
              thisKind.bonds[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == bondPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}
void PSFOutput::PrintAnglesInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxAngles[b], angleHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.angles.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d", thisKind.angles[i].a0 + atomID,
              thisKind.angles[i].a1 + atomID, thisKind.angles[i].a2 + atomID);
      ++lineEntry;
      if (lineEntry == anglePerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}
void PSFOutput::PrintDihedralsInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxDihs[b], dihedralHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.dihedrals.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d%8d", thisKind.dihedrals[i].a0 + atomID,
              thisKind.dihedrals[i].a1 + atomID,
              thisKind.dihedrals[i].a2 + atomID,
              thisKind.dihedrals[i].a3 + atomID);
      ++lineEntry;
      if (lineEntry == dihPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintImpropersInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxImps[b], improperHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.impropers.size(); ++i) {
      fprintf(outfile, "%8d%8d%8d%8d", thisKind.impropers[i].a0 + atomID,
              thisKind.impropers[i].a1 + atomID,
              thisKind.impropers[i].a2 + atomID,
              thisKind.impropers[i].a3 + atomID);
      ++lineEntry;
      if (lineEntry == impsPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintDonorsInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxDons[b], donorHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.donors.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.donors[i].a0 + atomID,
              thisKind.donors[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == donorPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintAcceptorsInBox(FILE *outfile, uint b) const {
  fprintf(outfile, headerFormat, boxAccs[b], acceptorHeader);
  uint atomID = 1;
  uint lineEntry = 0;
  for (MoleculeLookup::box_iterator thisMol = molLookRef.BoxBegin(b);
       thisMol != molLookRef.BoxEnd(b); thisMol++) {
    const MolKind &thisKind = molKinds[molecules->kIndex[*thisMol]];
    for (uint i = 0; i < thisKind.acceptors.size(); ++i) {
      fprintf(outfile, "%8d%8d", thisKind.acceptors[i].a0 + atomID,
              thisKind.acceptors[i].a1 + atomID);
      ++lineEntry;
      if (lineEntry == acceptorPerLine) {
        lineEntry = 0;
        fputc('\n', outfile);
      }
    }
    atomID += thisKind.atoms.size();
  }
  fputs("\n\n", outfile);
}

void PSFOutput::PrintNAMDCompliantSuffixInBox(FILE *outfile) const {
  // fprintf(outfile, headerFormat, 0, improperHeader);
  // fputs("\n\n", outfile);
  // fprintf(outfile, headerFormat, 0, donorHeader);
  // fputs("\n\n", outfile);
  // fprintf(outfile, headerFormat, 0, acceptorHeader);
  // fputs("\n\n", outfile);
  fprintf(outfile, headerFormat, 0, excludedHeader);
  fputs("\n\n", outfile);
  fprintf(outfile, headerFormat, 0, groupHeader);
}
