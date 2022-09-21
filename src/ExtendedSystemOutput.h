/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCD_OUTPUT_H
#define DCD_OUTPUT_H

#include <cstring>
#include <iostream>
#include <string> //to store lines of finished data.
#include <vector> //for molecule string storage.

#include "BasicTypes.h" //For uint
#include "Coordinates.h"
#include "DCDlib.h"
#include "MoleculeKind.h"
#include "Molecules.h"
#include "OutputAbstracts.h"
#include "PDBSetup.h" //For atoms class
#include "StaticVals.h"
#include "Velocity.h"
#include "Writer.h"

class System;
namespace config_setup {
struct Output;
}
class MoveSettings;
class MoleculeLookup;

struct ExtendedSystemOutput : OutputableBase {
public:
  ExtendedSystemOutput(System &sys, StaticVals const &statV);

  ~ExtendedSystemOutput() {
    if (x) {
      delete[] x; // y and z are continue of x
    }

    for (uint b = 0; b < BOX_TOTAL; ++b) {
      if (enableOut) {
        close_dcd_write(stateFileFileid[b]);
      }
      if (restartCoor[b]) {
        delete[] restartCoor[b];
      }
      if (restartVel[b]) {
        delete[] restartVel[b];
      }
      if (outDCDStateFile[b]) {
        delete[] outDCDStateFile[b];
      }
      if (outDCDRestartFile[b]) {
        delete[] outDCDRestartFile[b];
      }
      if (outVelRestartFile[b]) {
        delete[] outVelRestartFile[b];
      }
      if (outXSTFile[b]) {
        delete[] outXSTFile[b];
      }
      if (outXSCFile[b]) {
        delete[] outXSCFile[b];
      }
    }
  }

  // PDB does not need to sample on every step, so does nothing.
  virtual void Sample(const ulong step) {}

  virtual void Init(pdb_setup::Atoms const &atoms,
                    config_setup::Output const &output);

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);

private:
  // Copy cell length and angles to unitcell[6]
  void Copy_lattice_to_unitcell(double *unitcell, int box);
  // Unwrap and save coordinates of molecule in box, into *x, *y, *z
  void SetCoordinates(std::vector<int> &molInBox, const int box);
  // Unwrap and save coordinates of molecule in box, into *restartCoor
  // Save velocities of molecule in box, into *restartVel
  void SetMolInBox(const int box);
  // Return a vector that defines the box id for each molecule
  void SetMolBoxVec(std::vector<int> &mBox);
  // returns the total number of atoms in box
  int NumAtomInBox(const int box);
  // Write a binary restart file with coordinates
  void Write_binary_file(char *fname, int n, XYZ *vec);
  // Write header for dcd coordinate file
  void WriteDCDHeader(const int numAtoms, const int box);
  // Write header for xst and xsc files
  void Write_Extension_System_Header(Writer &outFile);
  // Write the cell basis info into xst and xsc files
  void Write_Extension_System_Data(Writer &outFile, const ulong step,
                                   const int box);

  MoveSettings &moveSetRef;
  MoleculeLookup &molLookupRef;
  BoxDimensions &boxDimRef;
  Molecules const &molRef;
  Coordinates &coordCurrRef;
  Velocity &velCurrRef;
  COM &comCurrRef;

  char *outDCDStateFile[BOX_TOTAL];
  char *outDCDRestartFile[BOX_TOTAL];
  char *outVelRestartFile[BOX_TOTAL];
  char *outXSTFile[BOX_TOTAL];
  char *outXSCFile[BOX_TOTAL];
  int stateFileFileid[BOX_TOTAL];
  bool outputVelocity; // output Velocity or not
  float *x, *y, *z;
  // AOS for restart binary format. NAMD internal data structure
  // is array of vector(XYZ)
  XYZ *restartCoor[BOX_TOTAL];
  XYZ *restartVel[BOX_TOTAL];
  // for extension system files
  Writer xstFile[BOX_TOTAL];
  Writer xscFile[BOX_TOTAL];
};

#endif /*DCD_OUTPUT_H*/
