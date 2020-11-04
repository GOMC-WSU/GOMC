/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCDOutput.h"              //For spec;
#include "EnsemblePreprocessor.h"   //For BOX_TOTAL, ensemble
#include "System.h"                 //for init
#include "StaticVals.h"             //for init
#include "MoleculeLookup.h"         //for lookup array (to get kind cnts, etc.)
#include "MoleculeKind.h"           //For kind names
#include "MoveSettings.h"           //For move settings/state
#include "PDBConst.h"               //For field locations/lengths
#include "StrStrmLib.h"             //For conversion from uint to string
#include <iostream>                 //for cout;

DCDOutput::DCDOutput(System  & sys, StaticVals const& statV) :
  moveSetRef(sys.moveSettings), molLookupRef(sys.molLookupRef),
  coordCurrRef(sys.coordinates), comCurrRef(sys.com),
  boxDimRef(sys.boxDimRef), molRef(statV.mol)
{
  x = NULL;
  y = NULL;
  z = NULL;
  enableStateOut = enableRestartOut = false;
  for(int b = 0; b < BOX_TOTAL; b++){
    stateFileFileid[b] = 0;
    restartCoor[b] = NULL;
    outDCDStateFile[b] = NULL;
    outDCDRestartFile[b] = NULL;
    outXSTFile[b] = NULL;
    outXSCFile[b] = NULL;
  }
}

void DCDOutput::Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output)
{
  enableStateOut = output.state_dcd.settings.enable;
  enableRestartOut = output.restart_dcd.settings.enable;
  enableOut = enableStateOut | enableRestartOut;
  stepsStatePerOut = output.state_dcd.settings.frequency;
  stepsRestartPerOut = output.restart_dcd.settings.frequency;
  if (stepsStatePerOut < stepsRestartPerOut) {
    stepsPerOut = stepsStatePerOut;
  } else {
    stepsPerOut = stepsRestartPerOut;
  }
  bool printNotify;
#ifndef NDEBUG
      printNotify = true;
#else
      printNotify = false;
#endif

  // Output dcd coordinates and xst file
  if (enableStateOut) {
    int numAtoms = coordCurrRef.Count();
    x = new float [3 * numAtoms];
    y = x + numAtoms;
    z = x + 2 * numAtoms;
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      std::string fileName = output.state_dcd.files.dcd.name[b];
      int baselen = strlen(fileName.c_str());
      outDCDStateFile[b] = new char[baselen];
      strcpy(outDCDStateFile[b], fileName.c_str());
      //  Write out the header with lattice parameter
      WriteDCDHeader(numAtoms, b);
      // prepare the xst file
      fileName = output.statistics.settings.uniqueStr.val;
      fileName += "_BOX_" + std::to_string(b) + ".xst";
      baselen = strlen(fileName.c_str());
      outXSTFile[b] = new char[baselen];
      strcpy(outXSTFile[b], fileName.c_str());
      xstFile[b].Init(fileName, " output XST", true, printNotify);
      xstFile[b].open();
      Write_Extention_System_Header(xstFile[b], b);
    }
  }

  // Output restart binary coordinates and xsc file
  if (enableRestartOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      std::string fileName = output.restart_dcd.files.dcd.name[b];
      restartCoor[b] = new XYZ[NumAtomInBox(b)];
      int baselen = strlen(fileName.c_str());
      outDCDRestartFile[b] = new char[baselen];
      strcpy(outDCDRestartFile[b], fileName.c_str());
      // prepare the xsc file
      fileName = output.statistics.settings.uniqueStr.val;
      fileName += "_BOX_" + std::to_string(b) + ".xsc";
      baselen = strlen(fileName.c_str());
      outXSCFile[b] = new char[baselen];
      strcpy(outXSCFile[b], fileName.c_str());
      xscFile[b].Init(fileName, " output XSC", true, printNotify);
    }
  }

  DoOutput(0);
}


void DCDOutput::Write_Extention_System_Header(Writer &outFile, const int box)
{
  outFile.file << "#$LABELS step";
  outFile.file << " a_x a_y a_z";
  outFile.file << " b_x b_y b_z";
  outFile.file << " c_x c_y c_z";
  outFile.file << " o_x o_y o_z";
  outFile.file << std::endl;
}


void DCDOutput::Write_Extention_System_data(Writer &outFile,
    const ulong step, const int box)
{
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

  outFile.file << " " << a.x;
  outFile.file << " " << a.y;
  outFile.file << " " << a.z;
  outFile.file << " " << b.x;
  outFile.file << " " << b.y;
  outFile.file << " " << b.z;
  outFile.file << " " << c.x;
  outFile.file << " " << c.y;
  outFile.file << " " << c.z;
  // Our origin is fix at origin
  outFile.file  << " " << 0 << " " << 0 << " " << 0;
  outFile.file << std::endl;
}

void DCDOutput::WriteDCDHeader(const int numAtoms, const int box)
{
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
  NSAVC = stepsStatePerOut;
  NPRIV = 0;
  NSTEP = NPRIV - NSAVC;
  NFILE = 0;
  int ret_code = write_dcdheader(stateFileFileid[box], 
      outDCDStateFile[box], numAtoms, NFILE, NPRIV, NSAVC, NSTEP,
      1.0/48.88841, 1);

  if (ret_code < 0) {
    char err_msg[257];
    sprintf(err_msg, "Writing of DCD header file %s failed!",
      outDCDStateFile[box]);
    NAMD_err(err_msg);
  }
}

void DCDOutput::DoOutput(const ulong step)
{
  // Output dcd coordinates and xst file
  if(enableStateOut) {
    int numAtoms = coordCurrRef.Count();
    // Determin which molecule is in which box. Assume we are in NVT
    // or NPT, otherwise, SetMolBoxVec would adjust the value.
    std::vector<int> molInBox(molRef.count, 0);
    SetMolBoxVec(molInBox);
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      //  Copy the coordinates for output
      SetCoordinates(x, y, z, molInBox, b);
      //  Write out the values for this step
      printf("Writing DCD coordinate to file %s at step %ld \n",
        outDCDStateFile[b], step+1);
      fflush(stdout);

      double unitcell[6];
      Copy_lattice_to_unitcell(unitcell, b);
      int ret_code = write_dcdstep(stateFileFileid[b], numAtoms, x, y, z, unitcell);
    
      if (ret_code < 0) {
        char err_msg[257];
        sprintf(err_msg, "Writing of DCD coordinate %s failed at step %ld!",
          outDCDStateFile[b], step+1);
        NAMD_err(err_msg);
      }
      printf("Finished writing DCD coordinate to file %s at step %ld \n",
        outDCDStateFile[b], step+1);

      // write the cellbasis data to xst file
      Write_Extention_System_data(xstFile[b], step, b);
    }
  }

  // Output restart binary coordinates and xsc file
  if (((step + 1) % stepsRestartPerOut == 0) && enableRestartOut) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      int numAtomInBox = NumAtomInBox(b);
      // Copy the coordinate data for each box into AOS
      SetMolInBox(b);
      printf("Writing binary restart coordinate to file %s at step %ld \n",
        outDCDRestartFile[b], step+1);
      //  Generate a binary restart file
      Write_binary_file(outDCDRestartFile[b], numAtomInBox, restartCoor[b]);
      printf("Finished writing binary restart coordinate to file %s at step %ld \n",
        outDCDRestartFile[b], step+1);
      
      // write XSC file
      NAMD_backup_file(outXSCFile[b], ".BAK");
      xscFile[b].openOverwrite();
      Write_Extention_System_Header(xscFile[b], b);
      // write the cellbasis data to xst file
      Write_Extention_System_data(xscFile[b], step, b);
      xscFile[b].close();
    }
  }

}

int DCDOutput::NumAtomInBox(const int box)
{
  int numAtoms = 0;
  int totKind = molLookupRef.GetNumKind();
  for (int k = 0; k < totKind; ++k) {
    numAtoms += (molLookupRef.NumKindInBox(k, box) * molRef.NumAtoms(k));
  }

  return numAtoms;
}

void DCDOutput::SetMolInBox(const int box)
{
  #if ENSEMBLE == GCMC || ENSEMBLE == GEMC      
  if(restartCoor[box]) {
    delete [] restartCoor[box];
    restartCoor[box] = new XYZ[NumAtomInBox(box)];
  }
  #endif

  uint i = 0, pStart = 0, pEnd = 0;
  XYZ ref, coor;
  MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(box),
                              end = molLookupRef.BoxEnd(box);
  while (m != end) {
    molRef.GetRangeStartStop(pStart, pEnd, *m);
    ref = comCurrRef.Get(*m);
    for (int p = pStart; p < pEnd; ++p) {
      coor = coordCurrRef.Get(p);
      boxDimRef.UnwrapPBC(coor, box, ref);
      
      restartCoor[box][i].x = coor.x;
      restartCoor[box][i].y = coor.y;
      restartCoor[box][i].z = coor.z;
      ++i;
    }
    ++m;
  }
}

void DCDOutput::Write_binary_file(char *fname, int n, XYZ *vec) 
{
  char errmsg[256];
  int fd;    //  File descriptor
  int32 n32 = n;
  NAMD_backup_file(fname, ".BAK");
  fd = NAMD_open(fname);

  sprintf(errmsg, "Error on write to binary restart coordinate file %s", fname);

  //  Write out the number of atoms and the vectors
  NAMD_write(fd, (char *) &n32, sizeof(int32), errmsg);
  NAMD_write(fd, (char *) vec, sizeof(XYZ)*n, errmsg);

  NAMD_close(fd, fname);

}

void DCDOutput::SetCoordinates(float *x, float *y, float *d, 
    std::vector<int> &molInBox, const int box)
{
  uint p, m, pStart = 0, pEnd = 0;
  int numMolecules = molRef.count;
  XYZ ref, coor;
  #if ENSEMBLE == NVT || ENSEMBLE == NPT
    //Loop through all molecules
    for (m = 0; m < numMolecules; ++m) {
      molRef.GetRangeStartStop(pStart, pEnd, m);
      ref = comCurrRef.Get(m);
      for (p = pStart; p < pEnd; ++p) {
        coor = coordCurrRef.Get(p);
        boxDimRef.UnwrapPBC(coor, box, ref);
        
        x[p] = coor.x;
        y[p] = coor.y;
        z[p] = coor.z;
      }
    }
  #else
    bool inThisBox; 
    //Loop through all molecules
    for (m = 0; m < numMolecules; ++m) {
      molRef.GetRangeStartStop(pStart, pEnd, m);
      ref = comCurrRef.Get(m);
      inThisBox = (molInBox[m] == box);
      for (p = pStart; p < pEnd; ++p) {
        if (inThisBox) {
          coor = coordCurrRef.Get(p);
          boxDimRef.UnwrapPBC(coor, box, ref);
        } else {
          coor.Reset();
        }
        
        x[p] = coor.x;
        y[p] = coor.y;
        z[p] = coor.z;
      }
    }
  #endif

}

void DCDOutput::Copy_lattice_to_unitcell(double *unitcell, int box) {
  unitcell[0] = boxDimRef.GetAxis(box).x;
  unitcell[2] = boxDimRef.GetAxis(box).y;
  unitcell[5] = boxDimRef.GetAxis(box).z;
  unitcell[1] = boxDimRef.cosAngle[box][2]; // gamma: cos AB
  unitcell[3] = boxDimRef.cosAngle[box][1]; // beta:  cos AC
  unitcell[4] = boxDimRef.cosAngle[box][0]; // alpha: cos BC
}


void DCDOutput::SetMolBoxVec(std::vector<int> & mBox)
{
  #if ENSEMBLE == GCMC || ENSEMBLE == GEMC
    for (int b = 0; b < BOX_TOTAL; ++b) {
      MoleculeLookup::box_iterator m = molLookupRef.BoxBegin(b),
                                  end = molLookupRef.BoxEnd(b);
      while (m != end) {
        mBox[*m] = b;
        ++m;
      }
    }
  #endif
}
