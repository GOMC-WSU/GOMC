#ifndef PSFOUTPUT_H
#define PSFOUTPUT_H

#include <vector>
#include <string>

#include "../lib/BasicTypes.h"
#include "MolSetup.h"

struct Molecules;

class PSFOutput
{
public:
    PSFOutput(const Molecules& molecules, mol_setup::MolMap& molMap, 
        const std::vector<std::string>& kindNames);
    
    //Output PSF file to filename using default remarks
    void PrintPSF(const std::string& filename) const;

    //Output PSF file to filename, specifying remarks
    void PrintPSF(const std::string& filename, const std::vector<std::string>& remarks) const;

private:
    const Molecules* molecules;
    std::vector<mol_setup::MolKind> molKinds;
    std::vector<std::string> molNames;
    uint totalAtoms;
    uint totalBonds;
    uint totalAngles;
    uint totalDihs;

    void PrintRemarks(FILE* outfile, const std::vector<std::string>& remarks) const;
    void PrintAtoms(FILE* outfile) const;
    void PrintBonds(FILE* outfile) const;
    void PrintAngles(FILE* outfile) const;
    void PrintDihedrals(FILE* outfile) const;

    void CountMolecules();
};

#endif