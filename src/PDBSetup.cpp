/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <vector>
#include <map> //for function lookup table.

#include "StrLib.h" //for string comparison wrapper
#include "PDBSetup.h" //Corresponding header to this body
#include "FixedWidthReader.h" //For fixed width reader
#include "ConfigSetup.h" //For restart info
#include "MoveConst.h"
#include <stdlib.h> //for exit
#include <string> // for to_string
#include "ExtendedSystem.h"


#if BOX_TOTAL == 1
const std::string PDBSetup::pdbAlias[] = {"system PDB coordinate file"};
#else
const std::string PDBSetup::pdbAlias[] = {"box 0 PDB coordinate file",
                                          "box 1 PDB coordinate file"
                                         };
#endif

namespace pdb_setup
{
void Remarks::SetRestart(config_setup::RestartSettings const& r )
{
  restart = r.enable;
  recalcTrajectory = r.recalcTrajectory;
  recalcTrajectoryBinary = r.recalcTrajectoryBinary;
  restartFromXSC = r.restartFromXSCFile;
  restartFromBinary = r.restartFromBinaryCoorFile;

  for(uint b = 0; b < BOX_TOTAL; b++) {
    if(recalcTrajectory)
      /* If the user provides a binary file when recalcTraj is true
        ignore PDB file except for first frame.  Even then, 
        the first frame can be overwritten by a binary coords. */
      reached[b] = recalcTrajectoryBinary;
    else
      reached[b] = true;
  }
}
void Remarks::Read(FixedWidthReader & pdb)
{
  using namespace pdb_entry::remark::field;
  using namespace pdb_entry::cryst1::field;

  if(restart) {
    //check if GOMC is taged and read the max dis, rot, vol value
    std::string varName;
    pdb.Get(varName, name::POS)
    .Get(disp[currBox], dis::POS)
    .Get(rotate[currBox], rot::POS)
    .Get(vol[currBox], vol::POS);

    CheckGOMC(varName);
    
  }
  if(recalcTrajectory && !recalcTrajectoryBinary) {
    std::string varName;
    pdb.Get(varName, name::POS)
    .Get(frameNumber[currBox], frameNum::POS)
    .Get(step[currBox], stepsNum::POS);

    if(frameNumber[currBox] == targetFrame[currBox])
      reached[currBox] = true;

    CheckGOMC(varName);
  }
}

void Remarks::CheckGOMC(std::string const& varName)
{
  using namespace pdb_entry::remark::field;
  if (!str::compare(varName, name::STR_GOMC)) {
    std::cerr << "ERROR: "
              << "GOMC file's identifying tag "
              << "\"REMARK     GOMC\" is missing"
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Remarks::Clear()
{
  frameSteps.clear();
}

void Cryst1::Read(FixedWidthReader & pdb)
{
  XYZ temp;
  using namespace pdb_entry::cryst1::field;
  hasVolume[currBox] = true;
  pdb.Get(temp.x, x::POS)
  .Get(temp.y, y::POS)
  .Get(temp.z, z::POS)
  .Get(cellAngle[currBox][0], ang_alpha::POS)
  .Get(cellAngle[currBox][1], ang_beta::POS)
  .Get(cellAngle[currBox][2], ang_gamma::POS);
  axis.Set(currBox, temp);
}

void Atoms::SetRestart(config_setup::RestartSettings const& r )
{
  restart = r.enable;
  restartFromBinary = r.restartFromBinaryCoorFile;
  recalcTrajectory = r.recalcTrajectory;
  //recalcTrajectory = r.recalcTrajectory && !r.recalcTrajectoryBinary;
  recalcTrajectoryBinary = r.recalcTrajectoryBinary;
}

void Atoms::Assign(std::string const& resName,
                   const char l_chain, const double l_x,
                   const double l_y, const double l_z,
                   const double l_beta)
{
  //box.push_back((bool)(restart?(uint)(l_occ):currBox));
  beta.push_back(l_beta);
  box.push_back(currBox);
  ++numAtomsInBox[currBox];
  resNames.push_back(resName);
  chainLetter.push_back(l_chain);

  // push the coordinates of atoms to x, y, and z
  x.push_back(l_x);
  y.push_back(l_y);
  z.push_back(l_z);

  count++;
}

void Atoms::Read(FixedWidthReader & file)
{
  using namespace pdb_entry::atom;
  char l_chain;
  uint resNum;
  std::string resName, atomName;
  double l_x, l_y, l_z, l_occ, l_beta;
  file.Get(atomName, field::alias::POS)
  .Get(resName, field::res_name::POS)
  .Get(resNum, field::res_num::POS)
  .Get(l_chain, field::chain::POS).Get(l_x, field::x::POS)
  .Get(l_y, field::y::POS).Get(l_z, field::z::POS)
  .Get(l_occ, field::occupancy::POS)
  .Get(l_beta, field::beta::POS);
  /* In the rare case you wanted 0-convention PDB Trajectories */
  if(recalcTrajectory && (uint)l_occ != currBox  && !recalcTrajectoryBinary) {
    return;
  }
  Assign(resName, l_chain, l_x, l_y, l_z, l_beta);
}

void Atoms::Clear()
{
  chainLetter.clear();
  x.clear();
  y.clear();
  z.clear();
  beta.clear();
  box.clear();
  resNames.clear();
  count = 0;
}

} //end namespace pdb_setup

void PDBSetup::Init(config_setup::RestartSettings const& restart,
                    config_setup::InFiles const& inFiles,
                    std::string const*const name, uint frameNum)
{
  // Clear the vectors for both atoms and remarks in case Init was called
  // more than once
  atoms.Clear();
  remarks.Clear();

  std::map<std::string, FWReadableBase *>::const_iterator dataKind;
  remarks.SetRestart(restart);
  atoms.SetRestart(restart);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    std::string varName = "";
    remarks.SetBox(b);
    remarks.SetFrameNumber(b, frameNum);
    cryst.SetBox(b);
    atoms.SetBox(b);
    std::string alias;
    if(remarks.recalcTrajectory) {
      sstrm::Converter toStr;
      std::string numStr = "";
      toStr << frameNum;
      toStr >> numStr;
      alias = pdbAlias[b] + " frame " + numStr;
    } else {
      alias = pdbAlias[b];
    }
    pdb[b].SetData(name[b], alias);

    // Open PDB only once and stay there
    // instead of re-opening it for every frame
    // refer to issue #131
    if(frameNum == 1)
      pdb[b].open();

    while (pdb[b].Read(varName, pdb_entry::label::POS)) {
      //If end of frame, and this is the frame we wanted,
      //end read on this file
      if (remarks.reached[b] && str::compare(varName, pdb_entry::end::STR)) {
        break;
      }

      //Call reader function if remarks were reached,
      // or it is a remark
      dataKind = dataKinds.find(varName);
      if (dataKind != dataKinds.end() &&
          (remarks.reached[b] ||
           str::compare(dataKind->first, pdb_entry::label::REMARK))) {
        dataKind->second->Read(pdb[b]);
      }
    }
    // If the recalcTrajectory is true and reached was still false
    // it means we couldn't find a remark and hence have to exit with error
    if(!remarks.reached[b] && remarks.recalcTrajectory && !remarks.recalcTrajectoryBinary) {
      std::cerr << "Error: Recalculate Trajectory is active..." << std::endl
                << ".. and couldn't find remark in PDB/DCD file!" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout.width(40);
    std::cout << std::left << "Finished reading: ";
    std::cout << "\t" << name[b] << std::endl;
  }

  if(remarks.recalcTrajectoryBinary)
    InitBinaryTrajectory(inFiles);
}

std::vector<ulong> PDBSetup::GetFrameSteps(std::string const*const name)
{
  std::map<std::string, FWReadableBase *>::const_iterator dataKind;
  remarks.SetBox(mv::BOX0);
  FixedWidthReader pdb(name[mv::BOX0], pdbAlias[mv::BOX0]);
  pdb.open();
  std::string varName;
  uint count = 0;
  while (pdb.Read(varName, pdb_entry::label::POS)) {
    if(varName == pdb_entry::label::REMARK) {
      dataKind = dataKinds.find(varName);
      dataKind->second->Read(pdb);
      remarks.frameSteps.push_back(remarks.step[mv::BOX0]);
      count++;
    }
  }
  return remarks.frameSteps;
}

bool PDBSetup::GetBinaryTrajectoryBoolean(){
  bool defined = false;
  for (uint b = 0; b < BOX_TOTAL; ++b){
    defined |= binTraj[b].defined;
  }
  return defined;
}

std::vector<ulong> PDBSetup::GetFrameStepsFromBinary(uint startStep, config_setup::InFiles const& inFiles){
  std::vector<ulong> frameSteps;
  for (uint b = 0; b < BOX_TOTAL; ++b){
    if(inFiles.binaryTrajectory.defined[b]){
      for (int i=0; i< binTraj[b].dcd->nsets; i++) {
        frameSteps.push_back(binTraj[b].dcd->istart + i*binTraj[b].dcd->nsavc);
      }
    }
  } 
  return frameSteps;
}

int PDBSetup::LoadBinaryTrajectoryStep(uint frameNum){
  int atomOffset = 0;
  for (int b = 0; b < BOX_TOTAL; ++b){
    if(binTraj[b].defined){
      int rc = DCDPlugin::read_next_timestep(binTraj[b].v, binTraj[b].natoms, &binTraj[b].timestep);
      if (rc) {
        fprintf(stderr, "error in read_next_timestep on frame %d\n", (int)frameNum);
        return 1;
      }
      /*NAMD trajectories use absent convention for non-occupying molecules in a box */
      if (binTraj[b].natoms != atoms.count){
        if (b == 0)
          atomOffset = 0;
        else 
          atomOffset = atoms.numAtomsInBox[0];
        for (int i = 0; i < binTraj[b].dcd->natoms; ++i) {
          atoms.x[atomOffset + i] = (double) binTraj[b].dcd->x[i];
          atoms.y[atomOffset + i] = (double) binTraj[b].dcd->y[i];
          atoms.z[atomOffset + i] = (double) binTraj[b].dcd->z[i];
        }
      } else {
        /* GOMC trajectories use [0.0, 0.0, 0.0]-convention for non-occupying molecules */
        for (int i = 0; i < binTraj[b].dcd->natoms; ++i){
          if (binTraj[b].dcd->x[i] != 0.0
                && binTraj[b].dcd->y[i] != 0.0
                  && binTraj[b].dcd->z[i] != 0.0){
            atoms.x[atomOffset + i] = (double) binTraj[b].dcd->x[i];
            atoms.y[atomOffset + i] = (double) binTraj[b].dcd->y[i];
            atoms.z[atomOffset + i] = (double) binTraj[b].dcd->z[i];
          }
        }
      }
    }
  }
  return 0;
}

int PDBSetup::InitBinaryTrajectory(config_setup::InFiles const& inFiles){
  for (int b = 0; b < BOX_TOTAL; ++b){
    if(inFiles.binaryTrajectory.defined[b]){
      binTraj[b].defined = true;
      int natoms = 0;
      binTraj[b].v = DCDPlugin::open_dcd_read(inFiles.binaryTrajectory.name[b].c_str(), "dcd", &natoms);
      if (!binTraj[b].v) {
        fprintf(stderr, "main) open_dcd_read failed for file %s\n", inFiles.binaryTrajectory.name[b].c_str());
        return 1;
      }
      binTraj[b].natoms = natoms;
      binTraj[b].dcd = (DCDPlugin::dcdhandle *)binTraj[b].v;
      binTraj[b].sizeMB = ((binTraj[b].dcd->natoms * 3.0) * binTraj[b].dcd->nsets * 4.0) / (1024.0 * 1024.0);
      binTraj[b].totalMB += binTraj[b].sizeMB; 
      printf("main) file: %s\n", inFiles.binaryTrajectory.name[b].c_str());
      printf("  %d atoms, %d frames, size: %6.1fMB\n", binTraj[b].dcd->natoms, binTraj[b].dcd->nsets, binTraj[b].sizeMB);

      
      binTraj[b].timestep.coords = (float *)malloc(3*sizeof(float)*binTraj[b].natoms);

    } else {
      binTraj[b].defined = false;
    }
  }
  return 0;
}
