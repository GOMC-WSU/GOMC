/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef EXTENDED_SYSTEM_H
#define EXTENDED_SYSTEM_H

#include <vector>

#include "BasicTypes.h"           //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "GeomLib.h"              // for PI and Dot product
#include "InputAbstracts.h"       //For FWReadableBase
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "PDBSetup.h"
#include "Velocity.h"
#include "XYZArray.h" //For box dimensions.
#include "XYZArray.h"

namespace config_setup {
struct RestartSettings;
struct InFiles;
struct Input;
} // namespace config_setup

namespace pdb_setup {
struct Cryst1;
}

class ExtendedSystem {
public:
  ExtendedSystem();
  ~ExtendedSystem(){};
  // Equality operator for unit testing
  bool operator==(const ExtendedSystem &other);
  void Init(PDBSetup &pdb, Velocity &vel, config_setup::Input &inputFiles,
            MoleculeLookup &molLookup, Molecules &mols);

private:
  // Reads the xsc file and store/calculate cellBasis data
  void ReadExtendedSystem(const char *filename, const int box);
  // Updates the cellBasis data in pdb data structure
  void UpdateCellBasis(PDBSetup &pdb, const int box);
  // Reads the binary coordinates and updates the X Y Z coordinates in pdb data
  // structure
  void UpdateCoordinate(PDBSetup &pdb, config_setup::Input &inputFiles,
                        MoleculeLookup &molLookup, Molecules &mols);
  void ReadCoordinate(PDBSetup &pdb, config_setup::Input &inputFiles,
                      MoleculeLookup &molLookup, Molecules &mols);
  void UpdateMinMaxAtoms(PDBSetup &pdb, config_setup::Input &inputFiles,
                         MoleculeLookup &molLookup, Molecules &mols);
  // Reads the binary velocities and updates the X Y Z velocity data structure
  void UpdateVelocity(Velocity &vel, config_setup::Input &inputFiles,
                      MoleculeLookup &molLookup, Molecules &mols);
  void ReadVelocity(PDBSetup &pdb, config_setup::Input &inputFiles,
                    MoleculeLookup &molLookup, Molecules &mols);
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

  // Offset array of molecule lookup array
  int boxMoleculeOffset[BOX_TOTAL + 1];

  // Stores the binary coordinates of both boxes
  std::vector<XYZ> binaryCoor;

  // Stores the binary velocities of both boxes
  std::vector<XYZ> binaryVeloc;
};

#endif /*EXTENDED_SYSTEM_H*/
