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
    int cmIndex = 0;
    UpdateCoordinate(pdb, inputFiles, molLookup, mols, cmIndex);    
  }
  // Read the binary velocity and update the buffer
  if(inputFiles.restart.restartFromBinaryVelFile) {
    int cmIndex = 0;
    for(int b = 0; b < BOX_TOTAL; b++) {
      if(inputFiles.files.binaryVelInput.defined[b]) {
        std::string fName = inputFiles.files.binaryVelInput.name[b];
        UpdateVelocity(pdb, vel, fName.c_str(), b, molLookup, mols, cmIndex);
      }
    }
  }
}

void ExtendedSystem::UpdateCoordinate(PDBSetup &pdb, config_setup::Input inputFiles, MoleculeLookup & molLookup,
                                      Molecules & mols, int & cmIndex)
{
    uint p, d, trajectoryI, dataI, placementStart, placementEnd, dataStart, dataEnd;
  // We must read restart PDB, which hold correct
  // number atom info in each Box
  int numAtoms = 0;
  int boxStart[BOX_TOTAL];
  boxStart[0] = 0;
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryCoorInput.defined[b]) {
      if (mols.restartFromCheckpoint){
        numAtoms += molLookup.restartedNumAtomsInBox[b];
        if (b == 1)
          boxStart[1] = molLookup.restartedNumAtomsInBox[0];
      } else {
        numAtoms += pdb.atoms.numAtomsInBox[b];
        if (b == 1)
          boxStart[1] = pdb.atoms.numAtomsInBox[0];
      }
    }
  }
  std::vector<XYZ> binaryCoor(numAtoms);
  for(int b = 0; b < BOX_TOTAL; b++) {
    if(inputFiles.files.binaryCoorInput.defined[b]) {
      std::string fName = inputFiles.files.binaryCoorInput.name[b];  
      if (mols.restartFromCheckpoint)
        read_binary_file(fName.c_str(), &binaryCoor[boxStart[b]], molLookup.restartedNumAtomsInBox[b]);
      else
        read_binary_file(fName.c_str(), &binaryCoor[boxStart[b]], pdb.atoms.numAtomsInBox[b]);
    }
  }
  
  //find the starting index
  for (int mol = 0; mol < molLookup.molLookupCount; mol++){
    trajectoryI = molLookup.restartMoleculeIndices[mol];
    dataI = mol;
    //Loop through particles in mol.
    mols.GetRestartOrderedRangeStartStop(placementStart, placementEnd, trajectoryI);
    mols.GetRangeStartStop(dataStart, dataEnd, dataI);
    for (p = placementStart, d = dataStart; p < placementEnd; ++p, ++d) {
      pdb.atoms.x[d] = binaryCoor[p].x;
      pdb.atoms.y[d] = binaryCoor[p].y;
      pdb.atoms.z[d] = binaryCoor[p].z;
    }
  }
}

void ExtendedSystem::UpdateVelocity(PDBSetup &pdb, Velocity &vel, const char *filename, const int box,
                                    MoleculeLookup & molLookup, Molecules & mols, int & cmIndex)
{
  // We must read restart PDB, which hold correct
  // number atom info in each Box
  int numAtoms = pdb.atoms.numAtomsInBox[box];
  int moleculeOffset = 0;
  XYZ *binaryVel;
  binaryVel = new XYZ[numAtoms];
  read_binary_file(filename, binaryVel, numAtoms);
  //find the starting index

  for(; cmIndex < (int) molLookup.molLookupCount; cmIndex++) {
    if(moleculeOffset >= numAtoms) break;
    int currentMolecule = molLookup.molLookup[cmIndex];
    int numberOfAtoms = mols.start[currentMolecule + 1] - mols.start[currentMolecule];
    int atomDestinationStart = mols.start[currentMolecule];

    for(int atom = 0; atom < numberOfAtoms; atom++) {
      vel.x[atomDestinationStart + atom] = binaryVel[moleculeOffset + atom].x;
      vel.y[atomDestinationStart + atom] = binaryVel[moleculeOffset + atom].y;
      vel.z[atomDestinationStart + atom] = binaryVel[moleculeOffset + atom].z;
    }

    moleculeOffset += numberOfAtoms;
  }

  delete [] binaryVel;
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


