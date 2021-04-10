/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EXTENDED_SYSTEM
#define EXTENDED_SYSTEM

#include "InputAbstracts.h" //For FWReadableBase
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "XYZArray.h" //For box dimensions.
#include "PDBSetup.h"
#include "Molecules.h"
#include "MoleculeLookup.h"
#include "GeomLib.h" // for PI and Dot product
#include <vector>

namespace config_setup
{
  struct RestartSettings;
  struct InFiles;
  struct Input;
}

namespace pdb_setup
{
  struct Cryst1;
}


class ExtendedSystem  {
  public:
  ExtendedSystem();
  ~ExtendedSystem() {};
  void Init(PDBSetup &pdb, config_setup::Input inputFiles, MoleculeLookup & molLookup, Molecules & mols);
  private:
    // Reads the xsc file and store/calculate cellBasis data
    void ReadExtendedSystem(const char *filename, const int box);
    // Updates the cellBasis data in pdb data structure
    void UpdateCellBasis(PDBSetup &pdb, const int box);
    // Reads the binary coordinates and updates the X Y Z coordinates in pdb data structure
    void UpdateCoordinate(PDBSetup &pdb, const char *filename, const int box, MoleculeLookup & molLookup, Molecules & mols, int & cmIndex);
    void UpdateCoordinateHybrid(PDBSetup &pdb, const char *filename, const int box);
    // the time steps in xsc file
    ulong firstStep;
    // Center of cell, but GOMC always uses 0 center
    XYZ center[BOX_TOTAL];
    // Cell basis vector
    XYZArray cellBasis[BOX_TOTAL];
    // Cell length
    XYZArray axis;
    // Cos values of alpha, beta, gamma
    double cosAngle[BOX_TOTAL][3];
    // Cell angle (alpha, beta, gamma)
    double cellAngle[BOX_TOTAL][3];
    // Check to see if xsc is defined
    bool hasCellBasis[BOX_TOTAL];
};


#endif /*EXTENDED_SYSTEM*/
