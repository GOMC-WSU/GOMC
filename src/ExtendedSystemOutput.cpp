/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "ExtendedSystemOutput.h" //For spec;

#include <iostream> //for cout;

#include "EnsemblePreprocessor.h" //For BOX_TOTAL, ensemble
#include "MoleculeKind.h"         //For kind names
#include "MoleculeLookup.h"       //for lookup array (to get kind cnts, etc.)
#include "MoveSettings.h"         //For move settings/state
#include "PDBConst.h"             //For field locations/lengths
#include "StaticVals.h"           //for init
#include "StrStrmLib.h"           //For conversion from uint to string
#include "System.h"               //for init

ExtendedSystemOutput::ExtendedSystemOutput(System &sys, StaticVals const &statV)
    : moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
      boxDimRef(sys.boxDimRef), molRef(statV.mol), velCurrRef(sys.vel),
      coordCurrRef(sys.coordinates), comCurrRef(sys.com) {
  x = NULL;
  y = NULL;
  z = NULL;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    stateFileFileid[b] = 0;
    restartCoor[b] = NULL;
    restartVel[b] = NULL;
    outDCDStateFile[b] = NULL;
    outDCDRestartFile[b] = NULL;
    outVelRestartFile[b] = NULL;
    outXSTFile[b] = NULL;
    outXSCFile[b] = NULL;
  }
}

void ExtendedSystemOutput::Init(pdb_setup::Atoms const &atoms,
                                config_setup::Output const &output) {
  enableOut = output.state_dcd.settings.enable;
  enableRestOut = output.restart.settings.enable;
  stepsPerOut = output.state_dcd.settings.frequency;
  stepsRestPerOut = output.restart.settings.frequency;
  outputVelocity = output.restart_vel.settings.enable;
  if (stepsPerOut >= stepsRestPerOut) {
    stepsPerOut = stepsRestPerOut;
  }
  bool printNotify;
#ifndef NDEBUG
  printNotify = true;
#else
  printNotify = false;
#endif

  // Output dcd coordinates and xst file
  if (enableOut) {
    int numAtoms = coordCurrRef.Count();
    x = new float[3 * numAtoms];
    y = x + numAtoms;
    z = x + 2 * numAtoms;
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      std::string fileName = output.state_dcd.files.dcd.name[b];
      int baselen = strlen(fileName.c_str());
      outDCDStateFile[b] = new char[baselen + 1];
      strcpy(outDCDStateFile[b], fileName.c_str());
      //  Write out the header with lattice parameter
      WriteDCDHeader(numAtoms, b);
      // prepare the xst file
      fileName = output.statistics.settings.uniqueStr.val;
      fileName += "_BOX_" + std::to_string(b) + ".xst";
      baselen = strlen(fileName.c_str());
      outXSTFile[b] = new char[baselen + 1];
      strcpy(outXSTFile[b], fileName.c_str());
      xstFile[b].Init(fileName, " output XST", true, printNotify);
      xstFile[b].open();
      Write_Extension_System_Header(xstFile[b]);
    }
    DoOutput(0);
  }

  // Output restart binary coordinates and xsc file
  if (enableRestOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      // prepare coor file
      std::string fileName = output.restart_dcd.files.dcd.name[b];
      restartCoor[b] = new XYZ[NumAtomInBox(b)];
      int baselen = strlen(fileName.c_str());
      outDCDRestartFile[b] = new char[baselen + 1];
      strcpy(outDCDRestartFile[b], fileName.c_str());

      // prepare vel file
      if (outputVelocity) {
        std::string fileName = output.restart_vel.files.dcd.name[b];
        restartVel[b] = new XYZ[NumAtomInBox(b)];
        baselen = strlen(fileName.c_str());
        outVelRestartFile[b] = new char[baselen + 1];
        strcpy(outVelRestartFile[b], fileName.c_str());
      }

      // prepare the xsc file
      fileName = output.statistics.settings.uniqueStr.val;
      fileName += "_BOX_" + std::to_string(b) + "_restart.xsc";
      baselen = strlen(fileName.c_str());
      outXSCFile[b] = new char[baselen + 1];
      strcpy(outXSCFile[b], fileName.c_str());
      xscFile[b].Init(fileName, " output XSC", true, printNotify);
    }
  }
}

void ExtendedSystemOutput::Write_Extension_System_Header(Writer &outFile) {
  outFile.file << "#$LABELS step";
  outFile.file << " a_x a_y a_z";
  outFile.file << " b_x b_y b_z";
  outFile.file << " c_x c_y c_z";
  outFile.file << " o_x o_y o_z";
  outFile.file << std::endl;
}

void ExtendedSystemOutput::Write_Extension_System_Data(Writer &outFile,
                                                       const ulong step,
                                                       const int box) {
  outFile.file.precision(12);
  if (step == 0) {
    outFile.file << 0;
  } else {
    outFile.file << step + 1;
  }

  // because we normalized the cell basis, we need to reverse it
  XYZ a = boxDimRef.cellBasis[box].Get(0) * boxDimRef.GetAxis(box).x;
  XYZ b = boxDimRef.cellBasis[box].Get(1) * boxDimRef.GetAxis(box).y;
  XYZ c = boxDimRef.cellBasis[box].Get(2) * boxDimRef.GetAxis(box).z;
  XYZ center(0.5 * a.Length(), 0.5 * b.Length(), 0.5 * c.Length());

  outFile.file << " " << a.x;
  outFile.file << " " << a.y;
  outFile.file << " " << a.z;
  outFile.file << " " << b.x;
  outFile.file << " " << b.y;
  outFile.file << " " << b.z;
  outFile.file << " " << c.x;
  outFile.file << " " << c.y;
  outFile.file << " " << c.z;
  // Our origin is fix at (0, 0 ,0 ), but we set it to half box length
  // to be compatible with NAMD
  outFile.file << " " << center.x << " " << center.y << " " << center.z;
  outFile.file << std::endl;
}

void ExtendedSystemOutput::WriteDCDHeader(const int numAtoms, const int box) {
  printf("Opening DCD coordinate file: %s \n", outDCDStateFile[box]);
  stateFileFileid[box] = open_dcd_write(outDCDStateFile[box]);

  if (stateFileFileid[box] == DCD_FILEEXISTS) {
    char err_msg[257];
    sprintf(err_msg, "DCD Coordinate file %s already exists!",
            outDCDStateFile[box]);
    NAMD_warn(err_msg);
  } else if (stateFileFileid[box] < 0) {
    char err_msg[257];
    sprintf(err_msg, "Couldn't open DCD coordinate file %s",
            outDCDStateFile[box]);
    NAMD_err(err_msg);
  }
  int NSAVC, NFILE, NPRIV, NSTEP;
  NSAVC = stepsPerOut;
  NPRIV = 0;
  NSTEP = NPRIV - NSAVC;
  NFILE = 0;
  int ret_code =
      write_dcdheader(stateFileFileid[box], outDCDStateFile[box], numAtoms,
                      NFILE, NPRIV, NSAVC, NSTEP, 1.0 / 48.88841, 1);

  if (ret_code < 0) {
    char err_msg[257];
    sprintf(err_msg, "Writing of DCD header file %s failed!",
            outDCDStateFile[box]);
    NAMD_err(err_msg);
  }
}

void ExtendedSystemOutput::DoOutput(const ulong step) {
  // Output dcd coordinates and xst file
  GOMC_EVENT_START(1, GomcProfileEvent::DCD_OUTPUT);
  int numAtoms = coordCurrRef.Count();
  // Determine which molecule is in which box. Assume we are in NVT
  // or NPT, otherwise, SetMolBoxVec would adjust the value.
  std::vector<int> molInBox(molRef.count, 0);
  if (BOX_TOTAL > 1)
    SetMolBoxVec(molInBox);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    //  Copy the coordinates for output
    SetCoordinates(molInBox, b);
    //  Write out the values for this step
    printf("Writing DCD coordinate to file %s at step %lu \n",
           outDCDStateFile[b], step + 1);
    fflush(stdout);

    double unitcell[6];
    Copy_lattice_to_unitcell(unitcell, b);
    int ret_code =
        write_dcdstep(stateFileFileid[b], numAtoms, x, y, z, unitcell);

    if (ret_code < 0) {
      char err_msg[257];
      sprintf(err_msg, "Writing of DCD coordinate %s failed at step %lu!",
              outDCDStateFile[b], step + 1);
      NAMD_err(err_msg);
    }
    printf("Finished writing DCD coordinate to file %s at step %lu \n",
           outDCDStateFile[b], step + 1);

    // write the cellbasis data to xst file
    Write_Extension_System_Data(xstFile[b], step, b);
    GOMC_EVENT_STOP(1, GomcProfileEvent::DCD_OUTPUT);
  }
}

void ExtendedSystemOutput::DoOutputRestart(const ulong step) {
  // Output restart binary coordinates and xsc file
  GOMC_EVENT_START(1, GomcProfileEvent::DCD_RESTART_OUTPUT);
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    int numAtomInBox = NumAtomInBox(b);
    // Copy the coordinate data for each box into AOS
    SetMolInBox(b);
    printf("Writing binary restart coordinate to file %s at step %lu \n",
           outDCDRestartFile[b], step + 1);
    //  Generate a binary restart file
    Write_binary_file(outDCDRestartFile[b], numAtomInBox, restartCoor[b]);
    printf(
        "Finished writing binary restart coordinate to file %s at step %lu \n",
        outDCDRestartFile[b], step + 1);

    // output restart velocities
    if (outputVelocity) {
      // Update the velocity in box
      velCurrRef.UpdateVelocityInBox(b);
      printf("Writing binary restart velocity to file %s at step %lu \n",
             outVelRestartFile[b], step + 1);
      //  Generate a binary restart velocity file
      Write_binary_file(outVelRestartFile[b], numAtomInBox, restartVel[b]);
      printf(
          "Finished writing binary restart velocity to file %s at step %lu \n",
          outVelRestartFile[b], step + 1);
    }
    // write XSC file
    NAMD_backup_file(outXSCFile[b], ".BAK");
    xscFile[b].openOverwrite();
    Write_Extension_System_Header(xscFile[b]);
    // write the cellbasis data to xst file
    Write_Extension_System_Data(xscFile[b], step, b);
    xscFile[b].close();
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::DCD_RESTART_OUTPUT);
}

int ExtendedSystemOutput::NumAtomInBox(const int box) {
  int numAtoms = 0;
  int totKind = molLookupRef.GetNumKind();
  for (int k = 0; k < totKind; ++k) {
    numAtoms += (molLookupRef.NumKindInBox(k, box) * molRef.NumAtoms(k));
  }

  return numAtoms;
}

void ExtendedSystemOutput::SetMolInBox(const int box) {
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
  if (restartCoor[box]) {
    delete[] restartCoor[box];
    restartCoor[box] = new XYZ[NumAtomInBox(box)];
  }
  if (restartVel[box]) {
    delete[] restartVel[box];
    restartVel[box] = new XYZ[NumAtomInBox(box)];
  }
#endif

  uint i = 0, pStart = 0, pEnd = 0;
  XYZ ref, coor;
  MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(box),
                               end = molLookupRef.BoxEnd(box);
  while (m != end) {
    molRef.GetRangeStartStop(pStart, pEnd, *m);
    ref = comCurrRef.Get(*m);
    for (uint p = pStart; p < pEnd; ++p) {
      coor = coordCurrRef.Get(p);
      boxDimRef.UnwrapPBC(coor, box, ref);

      restartCoor[box][i].x = coor.x;
      restartCoor[box][i].y = coor.y;
      restartCoor[box][i].z = coor.z;

      if (outputVelocity) {
        restartVel[box][i] = velCurrRef.Get(p);
      }

      ++i;
    }
    ++m;
  }
}

void ExtendedSystemOutput::Write_binary_file(char *fname, int n, XYZ *vec) {
  char errmsg[256];
  int fd; //  File descriptor
  int32 n32 = n;
  NAMD_backup_file(fname, ".BAK");
  fd = NAMD_open(fname);

  sprintf(errmsg, "Error on write to binary restart coordinate file %s", fname);

  //  Write out the number of atoms and the vectors
  NAMD_write(fd, (char *)&n32, sizeof(int32), errmsg);
  NAMD_write(fd, (char *)vec, sizeof(XYZ) * n, errmsg);

  NAMD_close(fd, fname);
}

void ExtendedSystemOutput::SetCoordinates(std::vector<int> &molInBox,
                                          const int box) {
  uint d, dataStart, dataEnd, dataI;
  int numMolecules = molRef.count;
  XYZ ref, coor;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
#if ENSEMBLE == NVT || ENSEMBLE == NPT
    // Loop through all molecules
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                 end = molLookupRef.BoxEnd(b);
    while (m != end) {
      dataI = *m;
      molRef.GetRangeStartStop(dataStart, dataEnd, dataI);
      ref = comCurrRef.Get(dataI);
      for (d = dataStart; d < dataEnd; ++d) {
        coor = coordCurrRef.Get(d);
        boxDimRef.UnwrapPBC(coor, box, ref);
        x[d] = coor.x;
        y[d] = coor.y;
        z[d] = coor.z;
      }
      ++m;
    }
  }
#else
    bool inThisBox;
    // Loop through all molecules
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                 end = molLookupRef.BoxEnd(b);
    while (m != end) {
      dataI = *m;
      molRef.GetRangeStartStop(dataStart, dataEnd, dataI);
      ref = comCurrRef.Get(dataI);
      inThisBox = (molInBox[dataI] == box);
      for (d = dataStart; d < dataEnd; ++d) {
        if (inThisBox) {
          coor = coordCurrRef.Get(d);
          boxDimRef.UnwrapPBC(coor, box, ref);
        } else {
          coor.Reset();
        }
        x[d] = coor.x;
        y[d] = coor.y;
        z[d] = coor.z;
      }
      ++m;
    }
  }
#endif
}

void ExtendedSystemOutput::Copy_lattice_to_unitcell(double *unitcell, int box) {
  unitcell[0] = boxDimRef.GetAxis(box).x;
  unitcell[2] = boxDimRef.GetAxis(box).y;
  unitcell[5] = boxDimRef.GetAxis(box).z;
  unitcell[1] = boxDimRef.cosAngle[box][2]; // gamma: cos AB
  unitcell[3] = boxDimRef.cosAngle[box][1]; // beta:  cos AC
  unitcell[4] = boxDimRef.cosAngle[box][0]; // alpha: cos BC
}

void ExtendedSystemOutput::SetMolBoxVec(std::vector<int> &mBox) {
  for (int b = 0; b < BOX_TOTAL; ++b) {
    MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                 end = molLookupRef.BoxEnd(b);
    while (m != end) {
      mBox[*m] = b;
      ++m;
    }
  }
}
