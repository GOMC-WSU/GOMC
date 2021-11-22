/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <vector>
#include "StrLib.h" //for string comparison wrapper
#include "ExtendedSystem.h" //Corresponding header to this body
#include "FixedWidthReader.h" //For fixed width reader
#include "ConfigSetup.h" //For restart info
#include "DCDlib.h" // for Error output
#include "Velocity.h" // for velocity data
#include <stdlib.h> //for exit
#include <string> // for to_string

ExtendedSystem::ExtendedSystem()
{
  firstStep = 0;
  axis.Init(BOX_TOTAL);
  for(int b = 0; b < BOX_TOTAL; b++) {
    cellBasis[b] = XYZArray(3);
    hasCellBasis[b] = false;
    center[b].Reset();
  }
}

void ExtendedSystem::Init(PDBSetup &pdb, Velocity &vel,  config_setup::Input inputFiles,
                          MoleculeLookup & molLookup, Molecules & mols)
{
  // Read the extended system file and update the cellBasis data
  if(inputFiles.restart.restartFromXSCFile) {
    for(int b = 0; b < BOX_TOTAL; b++) {
      if(inputFiles.files.xscInput.defined[b]) {
        std::string fName = inputFiles.files.xscInput.name[b];
        ReadExtendedSystem(fName.c_str(), b);
        UpdateCellBasis(pdb, b);
      }
    }
  }
  // Read the binary coordinate and update the PDB coordinate
  if(inputFiles.restart.restartFromBinaryCoorFile) {
    ReadCoordinate(pdb, inputFiles, molLookup, mols);
    UpdateCoordinate(pdb, molLookup, mols);    
    UpdateMinMaxAtoms(pdb, molLookup, mols);
  }
  // Read the binary velocity and update the buffer
  if(inputFiles.restart.restartFromBinaryVelFile) {
    ReadVelocity(pdb, inputFiles, molLookup, mols);
    UpdateVelocity(vel, molLookup, mols);
  }
}

void ExtendedSystem::UpdateCoordinate(PDBSetup &pdb, 
                                      MoleculeLookup & molLookup,
                                      Molecules & mols)
{
  uint p, d, trajectoryI, dataI, placementStart, placementEnd, dataStart, dataEnd;
  //find the starting index
  for (int mol = 0; mol < molLookup.molLookupCount; mol++){
    trajectoryI = molLookup.restartMoleculeIndices[mol];
    dataI = mol;
    //Loop through particles in mol.
    mols.GetRangeStartStop(placementStart, placementEnd, trajectoryI);
    mols.GetRestartOrderedRangeStartStop(dataStart, dataEnd, dataI);
    for (p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
      pdb.atoms.x[p] = binaryCoor[d].x;
      pdb.atoms.y[p] = binaryCoor[d].y;
      pdb.atoms.z[p] = binaryCoor[d].z;
    }
  }
}

void ExtendedSystem::ReadCoordinate(PDBSetup &pdb, config_setup::Input inputFiles, MoleculeLookup & molLookup,
                                     Molecules & mols){
    // We must read restart PDB, which hold correct
  // number atom info in each Box
  numAtoms = 0;
  boxStart[0] = 0;
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryCoorInput.defined[b]) {
      if (mols.restartFromCheckpoint){
        numAtomsInBox[b] = molLookup.restartedNumAtomsInBox[b];
      } else {
        numAtomsInBox[b] = pdb.atoms.numAtomsInBox[b];
      }
    }
  }

  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryCoorInput.defined[b]) {
      numAtoms += numAtomsInBox[b];
      if (b == 1)
        boxStart[1] = numAtomsInBox[0];
    }
  }

  binaryCoor.clear();
  binaryCoor.resize(numAtoms);
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryCoorInput.defined[b]) {
      std::string fName = inputFiles.files.binaryCoorInput.name[b];  
      if (mols.restartFromCheckpoint)
        read_binary_file(fName.c_str(), &binaryCoor[boxStart[b]], molLookup.restartedNumAtomsInBox[b]);
      else
        read_binary_file(fName.c_str(), &binaryCoor[boxStart[b]], pdb.atoms.numAtomsInBox[b]);
    }
  }
}

void ExtendedSystem::UpdateMinMaxAtoms(PDBSetup &pdb, 
                                      MoleculeLookup & molLookup,
                                      Molecules & mols){
  for (uint b = 0; b < BOX_TOTAL; b++) {
    int stRange, endRange;
    // -1 because we want to exclude the last array index  
    // Box 0
    // [0, numAtomsBox0)
    // Box 1
    // [numAtomsBox0, numAtomsBox0 + numAtomsBox1)

    // To prevent segfault
    if (numAtomsInBox[b] == 0)
      return;

    if (b == 0){
      stRange = 0;
      endRange = numAtomsInBox[0] - 1;
    } else if (b == 1) {
      stRange = numAtomsInBox[0];
      endRange = stRange + numAtomsInBox[1] - 1;
    } else {
      std::cout << "Error: Only Box 0 and Box 1 supported!" << std::endl;
      exit(EXIT_FAILURE);
    }

    XYZArray binaryCoorSOA(binaryCoor);

    pdb.atoms.min[b].x = *std::min_element(binaryCoorSOA.x + stRange, binaryCoorSOA.x + endRange);
    pdb.atoms.min[b].y = *std::min_element(binaryCoorSOA.y + stRange, binaryCoorSOA.y + endRange);
    pdb.atoms.min[b].z = *std::min_element(binaryCoorSOA.z + stRange, binaryCoorSOA.z + endRange);
    pdb.atoms.max[b].x = *std::max_element(binaryCoorSOA.x + stRange, binaryCoorSOA.x + endRange);
    pdb.atoms.max[b].y = *std::max_element(binaryCoorSOA.y + stRange, binaryCoorSOA.y + endRange);
    pdb.atoms.max[b].z = *std::max_element(binaryCoorSOA.z + stRange, binaryCoorSOA.z + endRange);
  }
}



void ExtendedSystem::ReadVelocity(PDBSetup &pdb, config_setup::Input inputFiles, MoleculeLookup & molLookup,
                                     Molecules & mols){
    // We must read restart PDB, which hold correct
  // number atom info in each Box
  // We define local variables here, for case that
  // Binary Velocity is provided but not binary coordinates.
  int numAtoms = 0;
  int boxStart[BOX_TOTAL];
  int numAtomsInBox[BOX_TOTAL];
  boxStart[0] = 0;
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryVelInput.defined[b]) {
      if (mols.restartFromCheckpoint){
        numAtomsInBox[b] = molLookup.restartedNumAtomsInBox[b];
      } else {
        numAtomsInBox[b] = pdb.atoms.numAtomsInBox[b];
      }
    }
  }

  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryVelInput.defined[b]) {
      numAtoms += numAtomsInBox[b];
      if (b == 1)
        boxStart[1] = numAtomsInBox[0];
    }
  }
  binaryVeloc.clear();
  binaryVeloc.resize(numAtoms);
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryVelInput.defined[b]) {
      std::string fName = inputFiles.files.binaryVelInput.name[b];  
      if (mols.restartFromCheckpoint)
        read_binary_file(fName.c_str(), &binaryVeloc[boxStart[b]], molLookup.restartedNumAtomsInBox[b]);
      else
        read_binary_file(fName.c_str(), &binaryVeloc[boxStart[b]], pdb.atoms.numAtomsInBox[b]);
    }
  }
}

void ExtendedSystem::UpdateVelocity(Velocity & vel, 
                                    MoleculeLookup & molLookup,
                                    Molecules & mols)
{
  uint p, d, trajectoryI, dataI, placementStart, placementEnd, dataStart, dataEnd;
  //find the starting index
  for (int mol = 0; mol < molLookup.molLookupCount; mol++){
    trajectoryI = molLookup.restartMoleculeIndices[mol];
    dataI = mol;
    //Loop through particles in mol.
    mols.GetRangeStartStop(placementStart, placementEnd, trajectoryI);
    mols.GetRestartOrderedRangeStartStop(dataStart, dataEnd, dataI);
    for (p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
      vel.x[p] = binaryVeloc[d].x;
      vel.y[p] = binaryVeloc[d].y;
      vel.z[p] = binaryVeloc[d].z;
    }
  }
}

void ExtendedSystem::UpdateCellBasis(PDBSetup &pdb, const int box)
{
  if (hasCellBasis[box]) {
    pdb.cryst.hasCellBasis[box] = true;
    // Important to set to false, so BoxDim reads the cellBasis vector
    // and not cell length and angle
    pdb.cryst.hasVolume[box] = false;
    cellBasis[box].CopyRange(pdb.cryst.cellBasis[box], 0, 0, 3);
    pdb.cryst.axis.Set(box, axis.Get(box));
    //pdb.remarks.step[box] = firstStep;
    for(int i = 0; i < 3; i++) {
      pdb.cryst.cellAngle[box][i] = cellAngle[box][i];
    }
  }
}

void ExtendedSystem::ReadExtendedSystem(const char *filename, const int box)
{
  char msg[257];
  sprintf(msg, "Info: Reading extended system file %s \n", filename);
  std::cout << msg << std::endl;
     
  std::ifstream xscFile(filename);
  if ( ! xscFile ){
    sprintf(msg, "Unable to open extended system file %s!\n", filename);
    NAMD_die(msg);
  }

  char labels[1024];
  do {
    if ( ! xscFile ){
      sprintf(msg, "Reading extended system file %s! \n", filename);
      NAMD_die(msg);
    }
    xscFile.getline(labels,1023);
  } while ( strncmp(labels,"#$LABELS ",9) );

  int a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z;
  a_x = a_y = a_z = b_x = b_y = b_z = c_x = c_y = c_z = -1;
  int o_x, o_y, o_z;
  o_x = o_y = o_z = -1;

  int pos = 0;
  char *l_i = labels + 8;
  while ( *l_i ) {
    if ( *l_i == ' ' ) { 
      ++l_i; 
      continue; 
    }
    char *l_i2;
    for ( l_i2 = l_i; *l_i2 && *l_i2 != ' '; ++l_i2 );
    if ( (l_i2 - l_i) == 3 && (l_i[1] == '_') ) {
      if (l_i[0] == 'a' && l_i[2] == 'x') a_x = pos;
      if (l_i[0] == 'a' && l_i[2] == 'y') a_y = pos;
      if (l_i[0] == 'a' && l_i[2] == 'z') a_z = pos;
      if (l_i[0] == 'b' && l_i[2] == 'x') b_x = pos;
      if (l_i[0] == 'b' && l_i[2] == 'y') b_y = pos;
      if (l_i[0] == 'b' && l_i[2] == 'z') b_z = pos;
      if (l_i[0] == 'c' && l_i[2] == 'x') c_x = pos;
      if (l_i[0] == 'c' && l_i[2] == 'y') c_y = pos;
      if (l_i[0] == 'c' && l_i[2] == 'z') c_z = pos;
      if (l_i[0] == 'o' && l_i[2] == 'x') o_x = pos;
      if (l_i[0] == 'o' && l_i[2] == 'y') o_y = pos;
      if (l_i[0] == 'o' && l_i[2] == 'z') o_z = pos;
    }
    ++pos;
    l_i = l_i2;
  }

  int numpos = pos;

  XYZ cell[3], origin;

  for ( pos = 0; pos < numpos; ++pos ) {
    double tmp;
    xscFile >> tmp;
    if ( ! xscFile ){
      sprintf(msg, "Reading BOX %d extended system file %s! \n",
              box+1, filename);
      NAMD_die(msg);
    }
    if ( pos == 0 ) firstStep = ulong(tmp);
    if ( pos == a_x ) cell[0].x = tmp;
    if ( pos == a_y ) cell[0].y = tmp;
    if ( pos == a_z ) cell[0].z = tmp;
    if ( pos == b_x ) cell[1].x = tmp;
    if ( pos == b_y ) cell[1].y = tmp;
    if ( pos == b_z ) cell[1].z = tmp;
    if ( pos == c_x ) cell[2].x = tmp;
    if ( pos == c_y ) cell[2].y = tmp;
    if ( pos == c_z ) cell[2].z = tmp;
    if ( pos == o_x ) origin.x = tmp;
    if ( pos == o_y ) origin.y = tmp;
    if ( pos == o_z ) origin.z = tmp;
  }


  sprintf(msg, "Info: Finished reading extended system file %s \n", filename);
  std::cout << msg << std::endl;
  // Store the cellBasis Vector, and calculate the cell angles
  hasCellBasis[box] = true;
  center[box] = origin;
  cellBasis[box].Set(0, cell[0]);
  cellBasis[box].Set(1, cell[1]);
  cellBasis[box].Set(2, cell[2]);
  axis.Set(box, cell[0].Length(), cell[1].Length(),
          cell[2].Length());
  //Find Cosine Angle of alpha, beta and gamma
  cosAngle[box][0] = geom::Dot(cellBasis[box].Get(1), cellBasis[box].Get(2)) /
                    (axis.Get(box).y * axis.Get(box).z);
  cosAngle[box][1] = geom::Dot(cellBasis[box].Get(0), cellBasis[box].Get(2)) /
                    (axis.Get(box).x * axis.Get(box).z);
  cosAngle[box][2] = geom::Dot(cellBasis[box].Get(0), cellBasis[box].Get(1)) /
                    (axis.Get(box).x * axis.Get(box).y);
  
  // Avoid numerical error
  for(int i = 0; i < 3; i++) {
    if(cosAngle[box][i] > 1.0) {
      cosAngle[box][i] = 1.0;
    }
    if(cosAngle[box][i] < -1.0) {
      cosAngle[box][i] = -1.0;
    }
    cellAngle[box][i] = float(acos(cosAngle[box][i]) * 180.0f * (float) M_1_PI);
  }
}


