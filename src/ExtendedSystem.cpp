/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "ExtendedSystem.h" //Corresponding header to this body

#include <stdlib.h> //for exit

#include <string> // for to_string
#include <vector>

#include "ConfigSetup.h"      //For restart info
#include "DCDlib.h"           // for Error output
#include "FixedWidthReader.h" //For fixed width reader
#include "StrLib.h"           //for string comparison wrapper
#include "Velocity.h"         // for velocity data

ExtendedSystem::ExtendedSystem() {
  firstStep = 0;
  axis.Init(BOX_TOTAL);
  for (int b = 0; b < BOX_TOTAL; b++) {
    cellBasis[b] = XYZArray(3);
    hasCellBasis[b] = false;
    center[b].Reset();
  }
}

// Equality operator for unit testing
bool ExtendedSystem::operator==(const ExtendedSystem &other) {
  bool result = true;
  result &= (firstStep == other.firstStep);
  result &= (axis == other.axis);
  // These are cleared after use, so unneccessary
  // result &= (binaryCoor == other.binaryCoor);
  // result &= (binaryVeloc == other.binaryVeloc);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    result &= (center[b] == other.center[b]);
    result &= (cellBasis[b] == other.cellBasis[b]);
    result &= (hasCellBasis[b] == other.hasCellBasis[b]);
    for (uint a = 0; a < 3; a++) {
      result &= (cosAngle[b][a] == other.cosAngle[b][a]);
      result &= (cellAngle[b][a] == other.cellAngle[b][a]);
    }
  }
  for (uint b = 0; b < BOX_TOTAL + 1; b++) {
    result &= (boxMoleculeOffset[b] == other.boxMoleculeOffset[b]);
  }
  return result;
}

void ExtendedSystem::Init(PDBSetup &pdb, Velocity &vel,
                          config_setup::Input &inputFiles,
                          MoleculeLookup &molLookup, Molecules &mols) {
  if (!inputFiles.restart.enable) {
    return;
  }
  // Read the extended system file and update the cellBasis data
  if (inputFiles.restart.restartFromXSCFile) {
    for (int b = 0; b < BOX_TOTAL; b++) {
      if (inputFiles.files.xscInput.defined[b]) {
        std::string fName = inputFiles.files.xscInput.name[b];
        ReadExtendedSystem(fName.c_str(), b);
        UpdateCellBasis(pdb, b);
      }
    }
  }

  boxMoleculeOffset[0] = 0;
  for (int box = 0; box < BOX_TOTAL; ++box) {
    boxMoleculeOffset[box + 1] =
        boxMoleculeOffset[box] + molLookup.NumInBox(box);
  }

  // Read the binary coordinate and update the PDB coordinate
  if (inputFiles.restart.restartFromBinaryCoorFile) {
    binaryCoor.clear();
    binaryCoor.resize(pdb.atoms.beta.size());
    ReadCoordinate(pdb, inputFiles, molLookup, mols);
    UpdateCoordinate(pdb, inputFiles, molLookup, mols);
    UpdateMinMaxAtoms(pdb, inputFiles, molLookup, mols);
    binaryCoor.clear();
  }
  // Read the binary velocity and update the buffer
  if (inputFiles.restart.restartFromBinaryVelFile) {
    binaryVeloc.clear();
    binaryVeloc.resize(pdb.atoms.beta.size());
    ReadVelocity(pdb, inputFiles, molLookup, mols);
    UpdateVelocity(vel, inputFiles, molLookup, mols);
    binaryVeloc.clear();
  }
}

void ExtendedSystem::UpdateCoordinate(PDBSetup &pdb,
                                      config_setup::Input &inputFiles,
                                      MoleculeLookup &molLookup,
                                      Molecules &mols) {
  uint p, d, trajectoryI, dataI, placementStart, placementEnd, dataStart,
      dataEnd;
  // find the starting index
  for (uint box = 0; box < BOX_TOTAL; ++box) {
    if (inputFiles.files.binaryCoorInput.defined[box]) {
      // find the starting index
      for (int mol = boxMoleculeOffset[box]; mol < boxMoleculeOffset[box + 1];
           ++mol) {
        dataI = mol;
        if (mols.restartFromCheckpoint) {
          trajectoryI = molLookup.molLookup[mol];
          mols.GetRestartOrderedRangeStartStop(dataStart, dataEnd, dataI);
        } else {
          trajectoryI = mol;
          mols.GetRangeStartStop(dataStart, dataEnd, dataI);
        }
        mols.GetRangeStartStop(placementStart, placementEnd, trajectoryI);
        // Loop through particles in mol.
        for (p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
          pdb.atoms.x[p] = binaryCoor[d].x;
          pdb.atoms.y[p] = binaryCoor[d].y;
          pdb.atoms.z[p] = binaryCoor[d].z;
        }
      }
    }
  }
}

void ExtendedSystem::ReadCoordinate(PDBSetup &pdb,
                                    config_setup::Input &inputFiles,
                                    MoleculeLookup &molLookup,
                                    Molecules &mols) {
  for (int b = 0; b < BOX_TOTAL; b++) {
    if (inputFiles.files.binaryCoorInput.defined[b]) {
      std::string fName = inputFiles.files.binaryCoorInput.name[b];
      read_binary_file(fName.c_str(), &binaryCoor[pdb.atoms.boxAtomOffset[b]],
                       pdb.atoms.numAtomsInBox[b]);
    }
  }
}

void ExtendedSystem::UpdateMinMaxAtoms(PDBSetup &pdb,
                                       config_setup::Input &inputFiles,
                                       MoleculeLookup &molLookup,
                                       Molecules &mols) {
  XYZArray binaryCoorSOA(binaryCoor);
  for (uint b = 0; b < BOX_TOTAL; b++) {
    if (inputFiles.files.binaryCoorInput.defined[b]) {
      // To prevent segfault
      if (pdb.atoms.numAtomsInBox[b] == 0)
        return;
      int stRange, endRange;
      stRange = pdb.atoms.boxAtomOffset[b];
      endRange = pdb.atoms.boxAtomOffset[b + 1];

      pdb.atoms.min[b].x = *std::min_element(binaryCoorSOA.x + stRange,
                                             binaryCoorSOA.x + endRange);
      pdb.atoms.min[b].y = *std::min_element(binaryCoorSOA.y + stRange,
                                             binaryCoorSOA.y + endRange);
      pdb.atoms.min[b].z = *std::min_element(binaryCoorSOA.z + stRange,
                                             binaryCoorSOA.z + endRange);
      pdb.atoms.max[b].x = *std::max_element(binaryCoorSOA.x + stRange,
                                             binaryCoorSOA.x + endRange);
      pdb.atoms.max[b].y = *std::max_element(binaryCoorSOA.y + stRange,
                                             binaryCoorSOA.y + endRange);
      pdb.atoms.max[b].z = *std::max_element(binaryCoorSOA.z + stRange,
                                             binaryCoorSOA.z + endRange);
    }
  }
}

void ExtendedSystem::ReadVelocity(PDBSetup &pdb,
                                  config_setup::Input &inputFiles,
                                  MoleculeLookup &molLookup, Molecules &mols) {
  for (int b = 0; b < BOX_TOTAL; b++) {
    if (inputFiles.files.binaryVelInput.defined[b]) {
      std::string fName = inputFiles.files.binaryVelInput.name[b];
      read_binary_file(fName.c_str(), &binaryVeloc[pdb.atoms.boxAtomOffset[b]],
                       pdb.atoms.numAtomsInBox[b]);
    }
  }
}

void ExtendedSystem::UpdateVelocity(Velocity &vel,
                                    config_setup::Input &inputFiles,
                                    MoleculeLookup &molLookup,
                                    Molecules &mols) {
  uint p, d, trajectoryI, dataI, placementStart, placementEnd, dataStart,
      dataEnd;
  for (uint box = 0; box < BOX_TOTAL; ++box) {
    if (inputFiles.files.binaryVelInput.defined[box]) {
      // find the starting index
      for (int mol = boxMoleculeOffset[box]; mol < boxMoleculeOffset[box + 1];
           ++mol) {
        dataI = mol;
        if (mols.restartFromCheckpoint) {
          trajectoryI = molLookup.molLookup[mol];
          mols.GetRestartOrderedRangeStartStop(dataStart, dataEnd, dataI);
        } else {
          trajectoryI = mol;
          mols.GetRangeStartStop(dataStart, dataEnd, dataI);
        }
        mols.GetRangeStartStop(placementStart, placementEnd, trajectoryI);
        // Loop through particles in mol.
        for (p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
          vel.x[p] = binaryVeloc[d].x;
          vel.y[p] = binaryVeloc[d].y;
          vel.z[p] = binaryVeloc[d].z;
        }
      }
    }
  }
}

void ExtendedSystem::UpdateCellBasis(PDBSetup &pdb, const int box) {
  pdb.cryst.hasCellBasis[box] = true;
  // Important to set to false, so BoxDim reads the cellBasis vector
  // and not cell length and angle
  pdb.cryst.hasVolume[box] = false;
  cellBasis[box].CopyRange(pdb.cryst.cellBasis[box], 0, 0, 3);
  pdb.cryst.axis.Set(box, axis.Get(box));
  // pdb.remarks.step[box] = firstStep;
  for (int i = 0; i < 3; i++) {
    pdb.cryst.cellAngle[box][i] = cellAngle[box][i];
  }
}

void ExtendedSystem::ReadExtendedSystem(const char *filename, const int box) {
  char msg[257];
  sprintf(msg, "Info: Reading extended system file %s \n", filename);
  std::cout << msg << std::endl;

  std::ifstream xscFile(filename);
  if (!xscFile) {
    sprintf(msg, "Unable to open extended system file %s!\n", filename);
    NAMD_die(msg);
  }

  char labels[1024];
  do {
    if (!xscFile) {
      sprintf(msg, "Reading extended system file %s! \n", filename);
      NAMD_die(msg);
    }
    xscFile.getline(labels, 1023);
  } while (strncmp(labels, "#$LABELS ", 9));

  int a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z;
  a_x = a_y = a_z = b_x = b_y = b_z = c_x = c_y = c_z = -1;
  int o_x, o_y, o_z;
  o_x = o_y = o_z = -1;

  int pos = 0;
  char *l_i = labels + 8;
  while (*l_i) {
    if (*l_i == ' ') {
      ++l_i;
      continue;
    }
    char *l_i2;
    for (l_i2 = l_i; *l_i2 && *l_i2 != ' '; ++l_i2)
      ;
    if ((l_i2 - l_i) == 3 && (l_i[1] == '_')) {
      if (l_i[0] == 'a' && l_i[2] == 'x')
        a_x = pos;
      if (l_i[0] == 'a' && l_i[2] == 'y')
        a_y = pos;
      if (l_i[0] == 'a' && l_i[2] == 'z')
        a_z = pos;
      if (l_i[0] == 'b' && l_i[2] == 'x')
        b_x = pos;
      if (l_i[0] == 'b' && l_i[2] == 'y')
        b_y = pos;
      if (l_i[0] == 'b' && l_i[2] == 'z')
        b_z = pos;
      if (l_i[0] == 'c' && l_i[2] == 'x')
        c_x = pos;
      if (l_i[0] == 'c' && l_i[2] == 'y')
        c_y = pos;
      if (l_i[0] == 'c' && l_i[2] == 'z')
        c_z = pos;
      if (l_i[0] == 'o' && l_i[2] == 'x')
        o_x = pos;
      if (l_i[0] == 'o' && l_i[2] == 'y')
        o_y = pos;
      if (l_i[0] == 'o' && l_i[2] == 'z')
        o_z = pos;
    }
    ++pos;
    l_i = l_i2;
  }

  int numpos = pos;

  XYZ cell[3], origin;

  for (pos = 0; pos < numpos; ++pos) {
    double tmp;
    xscFile >> tmp;
    if (!xscFile) {
      sprintf(msg, "Reading BOX %d extended system file %s! \n", box + 1,
              filename);
      NAMD_die(msg);
    }
    if (pos == 0)
      firstStep = ulong(tmp);
    if (pos == a_x)
      cell[0].x = tmp;
    if (pos == a_y)
      cell[0].y = tmp;
    if (pos == a_z)
      cell[0].z = tmp;
    if (pos == b_x)
      cell[1].x = tmp;
    if (pos == b_y)
      cell[1].y = tmp;
    if (pos == b_z)
      cell[1].z = tmp;
    if (pos == c_x)
      cell[2].x = tmp;
    if (pos == c_y)
      cell[2].y = tmp;
    if (pos == c_z)
      cell[2].z = tmp;
    if (pos == o_x)
      origin.x = tmp;
    if (pos == o_y)
      origin.y = tmp;
    if (pos == o_z)
      origin.z = tmp;
  }

  sprintf(msg, "Info: Finished reading extended system file %s \n", filename);
  std::cout << msg << std::endl;
  // Store the cellBasis Vector, and calculate the cell angles
  hasCellBasis[box] = true;
  center[box] = origin;
  cellBasis[box].Set(0, cell[0]);
  cellBasis[box].Set(1, cell[1]);
  cellBasis[box].Set(2, cell[2]);
  axis.Set(box, cell[0].Length(), cell[1].Length(), cell[2].Length());
  // Find Cosine Angle of alpha, beta and gamma
  cosAngle[box][0] = geom::Dot(cellBasis[box].Get(1), cellBasis[box].Get(2)) /
                     (axis.Get(box).y * axis.Get(box).z);
  cosAngle[box][1] = geom::Dot(cellBasis[box].Get(0), cellBasis[box].Get(2)) /
                     (axis.Get(box).x * axis.Get(box).z);
  cosAngle[box][2] = geom::Dot(cellBasis[box].Get(0), cellBasis[box].Get(1)) /
                     (axis.Get(box).x * axis.Get(box).y);

  // Avoid numerical error
  for (int i = 0; i < 3; i++) {
    if (cosAngle[box][i] > 1.0) {
      cosAngle[box][i] = 1.0;
    }
    if (cosAngle[box][i] < -1.0) {
      cosAngle[box][i] = -1.0;
    }
    cellAngle[box][i] = float(acos(cosAngle[box][i]) * 180.0f * (float)M_1_PI);
  }
}
