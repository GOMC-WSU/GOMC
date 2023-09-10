/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "ConfigSetup.h"

#include <iomanip>
#include <map>    //for function handle storage.
#include <string> //for var names, etc.
#include <vector>

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623158e+308
#endif

ConfigSetup::ConfigSetup(void) {
  int i;
  exptMode = false;
  in.restart.enable = false;
  in.restart.step = ULONG_MAX;
  in.restart.recalcTrajectory = false;
  in.restart.restartFromCheckpoint = false;
  in.restart.restartFromBinaryCoorFile = false;
  in.restart.restartFromBinaryVelFile = false;
  in.restart.restartFromXSCFile = false;
  in.prng.seed = UINT_MAX;
  in.prngParallelTempering.seed = UINT_MAX;
  sys.elect.readEwald = false;
  sys.elect.readElect = false;
  sys.elect.readCache = false;
  sys.elect.ewald = false;
  sys.elect.enable = false;
  sys.elect.cache = false;
  sys.elect.tolerance = DBL_MAX;
  sys.elect.oneFourScale = DBL_MAX;
  sys.elect.dielectric = DBL_MAX;
  sys.memcVal.enable = false;
  sys.neMTMCVal.enable = false;
  sys.intraMemcVal.enable = false;
  sys.step.start = 0;
  sys.step.total = ULONG_MAX;
  sys.step.equil = ULONG_MAX;
  sys.step.adjustment = ULONG_MAX;
  sys.step.initStepRead = false;
  sys.step.initStep = ULONG_MAX;
  sys.step.pressureCalcFreq = ULONG_MAX;
  sys.step.pressureCalc = true;
  sys.step.parallelTempFreq = ULONG_MAX;
  sys.step.parallelTemperingAttemptsPerExchange = 0;
  sys.step.pressureCalc = false;
  in.ffKind.numOfKinds = 0;
  sys.exclude.EXCLUDE_KIND = UINT_MAX;
  in.prng.kind = "";
  for (i = 0; i < BOX_TOTAL; i++) {
    in.files.pdb.name[i] = "";
    in.files.psf.name[i] = "";
    in.files.binaryCoorInput.name[i] = "";
    in.files.binaryVelInput.name[i] = "";
    in.files.xscInput.name[i] = "";
  }
#if ENSEMBLE == GEMC
  sys.gemc.kind = UINT_MAX;
  sys.gemc.pressure = DBL_MAX;
#endif
#if ENSEMBLE == NPT
  sys.gemc.kind = mv::GEMC_NPT;
  sys.gemc.pressure = DBL_MAX;
#endif

  sys.T.inKelvin = DBL_MAX;
  sys.ff.VDW_KIND = UINT_MAX;
  sys.ff.doTailCorr = true;
  sys.ff.doImpulsePressureCorr = false;
  sys.ff.rswitch = DBL_MAX;
  sys.ff.cutoff = DBL_MAX;
  sys.ff.cutoffLow = DBL_MAX;
  sys.ff.vdwGeometricSigma = false;
  sys.moves.displace = DBL_MAX;
  sys.moves.rotate = DBL_MAX;
  sys.moves.intraSwap = DBL_MAX;
  sys.moves.multiParticleEnabled = false;
  sys.moves.multiParticle = DBL_MAX;
  sys.moves.multiParticleBrownian = DBL_MAX;
  sys.moves.regrowth = DBL_MAX;
  sys.moves.crankShaft = DBL_MAX;
  sys.moves.intraTargetedSwap = DBL_MAX;
  sys.moves.intraMemc = DBL_MAX;
  out.state.settings.enable = false;
  out.restart.settings.enable = false;
  out.state_dcd.settings.enable = false;
  out.restart_vel.settings.enable = false;
  out.console.enable = true;
  out.statistics.settings.block.enable = true;
#if ENSEMBLE == GCMC
  sys.chemPot.isFugacity = false;
  out.statistics.settings.hist.enable = false;
  out.statistics.settings.hist.frequency = ULONG_MAX;
  out.state.files.hist.histName = "";
  out.state.files.hist.letter = "";
  out.state.files.hist.number = "";
  out.state.files.hist.sampleName = "";
  out.state.files.hist.stepsPerHistSample = UINT_MAX;
#endif
  out.checkpoint.enable = false;
  out.checkpoint.frequency = ULONG_MAX;
  out.statistics.settings.uniqueStr.val = "";
  out.state.settings.frequency = ULONG_MAX;
  out.restart.settings.frequency = ULONG_MAX;
  out.state_dcd.settings.frequency = ULONG_MAX;
  out.console.frequency = ULONG_MAX;
  out.statistics.settings.block.frequency = ULONG_MAX;
  out.statistics.vars.energy.block = false;
  out.statistics.vars.energy.fluct = false;
  out.statistics.vars.pressure.block = false;
  out.statistics.vars.pressure.fluct = false;
  out.statistics.vars.surfaceTension.block = false;
  out.statistics.vars.surfaceTension.fluct = false;
#ifdef VARIABLE_PARTICLE_NUMBER
  sys.moves.transfer = DBL_MAX;
  sys.moves.memc = DBL_MAX;
  sys.moves.neMolTransfer = DBL_MAX;
  sys.moves.targetedSwap = DBL_MAX;
  sys.cbmcTrials.bonded.ang = UINT_MAX;
  sys.cbmcTrials.bonded.dih = UINT_MAX;
  sys.cbmcTrials.nonbonded.first = UINT_MAX;
  sys.cbmcTrials.nonbonded.nth = UINT_MAX;
  out.statistics.vars.molNum.block = false;
  out.statistics.vars.molNum.fluct = false;
  sys.volume.cstVolBox0 = false;
#endif
#ifdef VARIABLE_VOLUME
  sys.moves.volume = DBL_MAX;
  out.statistics.vars.volume.block = false;
  out.statistics.vars.volume.fluct = false;
#endif
  out.statistics.vars.density.block = false;
  out.statistics.vars.density.fluct = false;
}

int ConfigSetup::stringtoi(const std::string &s) {
  std::istringstream str(s);
  uint i;
  str >> i;
  return i;
}

double ConfigSetup::stringtod(const std::string &s) {
  std::istringstream str(s);
  double i;
  str >> i;
  return i;
}

bool ConfigSetup::checkBool(std::string str) {
  // short circuit for long strings
  if (str.length() > 5) {
    std::cout << "ERROR: " << str << "couldn't be recognized!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capitalize string
  for (uint k = 0; k < str.length(); k++) {
    str[k] = toupper(str[k]);
  }

  if (str == "ON" || str == "TRUE" || str == "YES")
    return true;
  else if (str == "OFF" || str == "FALSE" || str == "NO")
    return false;
  std::cout << "ERROR: " << str << "couldn't be recognized!" << std::endl;
  exit(EXIT_FAILURE);
}

// Same as checkBool but doesn't exit with an ERROR if the string doesn't
// represent a boolean value. Use for cases where we aren't sure the argument
// string should be a boolean.
bool ConfigSetup::isBool(std::string str) {
  // short circuit for long strings
  if (str.length() > 5)
    return false;

  // capitalize string
  for (uint k = 0; k < str.length(); k++) {
    str[k] = toupper(str[k]);
  }

  return (str == "ON" || str == "TRUE" || str == "YES" || str == "OFF" ||
          str == "FALSE" || str == "NO");
}

bool ConfigSetup::CheckString(std::string str1, std::string str2) {
  for (uint k = 0; k < str1.length(); k++) {
    str1[k] = toupper(str1[k]);
  }

  for (uint j = 0; j < str2.length(); j++) {
    str2[j] = toupper(str2[j]);
  }

  return (str1 == str2);
}

void ConfigSetup::Init(const char *fileName, MultiSim const *const &multisim) {
  std::vector<std::string> line;

  reader.Open(fileName);
  printf("\n%-40s %-s\n", "Reading Input File:", fileName);
  while (reader.readNextLine(line)) {
    if (line.size() == 0)
      continue;

    if (CheckString(line[0], "Restart")) {
      in.restart.enable = checkBool(line[1]);
      if (in.restart.enable) {
        printf("%-40s %-s \n", "Info: Restart simulation", "Active");
      }
    } else if (CheckString(line[0], "FirstStep")) {
      in.restart.step = stringtoi(line[1]);
    } else if (CheckString(line[0], "PRNG")) {
      if (line[1] == "RANDOM") {
        printf("%-40s %-s \n", "Info: Random seed", "Active");
        in.prng.kind = line[1];
      } else if (line[1] == "INTSEED") {
        printf("%-40s %-s \n", "Info: Integer seed", "Active");
        in.prng.kind = line[1];
      } else {
        std::cout << "ERROR: PRNG can only be \"RANDOM\" or \"INTSEED\"!\n";
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "PRNG_ParallelTempering")) {
      in.prngParallelTempering.kind = line[1];
      if ("RANDOM" == line[1])
        printf("%-40s %-s \n", "Info: Random seed", "Active");
    } else if (CheckString(line[0], "ParaTypeCHARMM")) {
      if (checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isCHARMM = true;
        printf("%-40s %-s \n", "Info: PARAMETER file", "CHARMM format");
      }
    } else if (CheckString(line[0], "ParaTypeEXOTIC")) {
      if (checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isCHARMM = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isEXOTIC = true;
        printf("%-40s %-s \n", "Info: PARAMETER file", "MIE format");
      }
    } else if (CheckString(line[0], "ParaTypeMIE")) {
      if (checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isCHARMM = false;
        in.ffKind.isMARTINI = false;
        in.ffKind.isEXOTIC = true;
        printf("%-40s %-s \n", "Info: PARAMETER file", "MIE format");
      }
    } else if (CheckString(line[0], "ParaTypeMARTINI")) {
      if (checkBool(line[1])) {
        in.ffKind.numOfKinds++;
        in.ffKind.isEXOTIC = false;
        in.ffKind.isMARTINI = true;
        in.ffKind.isCHARMM = true;
        printf("%-40s %-s \n", "Info: PARAMETER file",
               "MARTINI using CHARMM format");
      }
    } else if (CheckString(line[0], "Parameters")) {
      if (line.size() > 1) {
        in.files.param.push_back(config_setup::FileName(line[1], true));
      }
    } else if (CheckString(line[0], "Coordinates")) {
      uint boxnum = stringtoi(line[1]);
      if (boxnum >= BOX_TOTAL) {
        std::cout << "ERROR: Simulation requires " << BOX_TOTAL
                  << " PDB file(s)!\n";
        exit(EXIT_FAILURE);
      }
      if (multisim != NULL) {
        in.files.pdb.name[boxnum] =
            multisim->replicaInputDirectoryPath + line[2];
      } else {
        in.files.pdb.name[boxnum] = line[2];
      }
      in.files.pdb.defined[boxnum] = true;
    } else if (CheckString(line[0], "Structure")) {
      uint boxnum = stringtoi(line[1]);
      if (boxnum >= BOX_TOTAL) {
        std::cout << "ERROR: Simulation requires " << BOX_TOTAL
                  << " PSF file(s)!\n";
        exit(EXIT_FAILURE);
      }
      if (multisim != NULL) {
        in.files.psf.name[boxnum] =
            multisim->replicaInputDirectoryPath + line[2];
      } else {
        in.files.psf.name[boxnum] = line[2];
      }
      in.files.psf.defined[boxnum] = true;
    } else if (CheckString(line[0], "binCoordinates")) {
      uint boxnum = stringtoi(line[1]);
      if (boxnum >= BOX_TOTAL) {
        std::cout << "ERROR: Simulation requires " << BOX_TOTAL
                  << " binary coordinate file(s)!\n";
        exit(EXIT_FAILURE);
      }
      if (multisim != NULL) {
        in.files.binaryCoorInput.name[boxnum] =
            multisim->replicaInputDirectoryPath + line[2];
      } else {
        in.files.binaryCoorInput.name[boxnum] = line[2];
      }
      in.files.binaryCoorInput.defined[boxnum] = true;
      in.restart.restartFromBinaryCoorFile = true;
    } else if (CheckString(line[0], "binVelocities")) {
      uint boxnum = stringtoi(line[1]);
      if (boxnum >= BOX_TOTAL) {
        std::cout << "ERROR: Simulation requires " << BOX_TOTAL
                  << " binary velocity file(s)!\n";
        exit(EXIT_FAILURE);
      }
      if (multisim != NULL) {
        in.files.binaryVelInput.name[boxnum] =
            multisim->replicaInputDirectoryPath + line[2];
      } else {
        in.files.binaryVelInput.name[boxnum] = line[2];
      }
      in.files.binaryVelInput.defined[boxnum] = true;
      in.restart.restartFromBinaryVelFile = true;
      // If we read the binary vel, we also output the restart vel. Otherwise
      // we dont have any output velocity
      out.restart_vel.settings.enable = true;
    } else if (CheckString(line[0], "extendedSystem")) {
      uint boxnum = stringtoi(line[1]);
      if (boxnum >= BOX_TOTAL) {
        std::cout << "ERROR: Simulation requires " << BOX_TOTAL
                  << " extended system file(s)!\n";
        exit(EXIT_FAILURE);
      }
      if (multisim != NULL) {
        in.files.xscInput.name[boxnum] =
            multisim->replicaInputDirectoryPath + line[2];
      } else {
        in.files.xscInput.name[boxnum] = line[2];
      }
      in.files.xscInput.defined[boxnum] = true;
      in.restart.restartFromXSCFile = true;
    } else if (CheckString(line[0], "Checkpoint")) {
      if (line.size() == 3) {
        in.files.checkpoint.defined[0] = in.restart.restartFromCheckpoint =
            checkBool(line[1]);
        if (multisim != NULL) {
          in.files.checkpoint.name[0] =
              multisim->replicaInputDirectoryPath + line[2];
        } else {
          in.files.checkpoint.name[0] = line[2];
        }
      } else {
        if (isBool(line[1])) {
          in.files.checkpoint.defined[0] = in.restart.restartFromCheckpoint =
              checkBool(line[1]);
          in.files.checkpoint.name[0] = "";
        } else {
          in.files.checkpoint.defined[0] = in.restart.restartFromCheckpoint =
              true;
          if (multisim != NULL) {
            in.files.checkpoint.name[0] =
                multisim->replicaInputDirectoryPath + line[1];
          } else {
            in.files.checkpoint.name[0] = line[1];
          }
        }
      }
    }
#if ENSEMBLE == GEMC
    else if (CheckString(line[0], "GEMC")) {
      if (CheckString(line[1], "NVT")) {
        sys.gemc.kind = mv::GEMC_NVT;
        printf("Info: Running NVT_GEMC\n");
      } else if (CheckString(line[1], "NPT")) {
        sys.gemc.kind = mv::GEMC_NPT;
        printf("Info: Running NPT_GEMC\n");
      }
    } else if (CheckString(line[0], "Pressure")) {
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-40s %-4.4f bar\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
#if ENSEMBLE == NPT
    else if (CheckString(line[0], "Pressure")) {
      sys.gemc.kind = mv::GEMC_NPT;
      sys.gemc.pressure = stringtod(line[1]);
      printf("%-40s %-4.4f bar\n", "Info: Input Pressure", sys.gemc.pressure);
      sys.gemc.pressure *= unit::BAR_TO_K_MOLECULE_PER_A3;
    }
#endif
    else if (CheckString(line[0], "Temperature")) {
      if (line.size() > 2 && multisim != NULL) {
        sys.T.inKelvin = stringtod(line[multisim->worldRank + 1]);
      } else {
        sys.T.inKelvin = stringtod(line[1]);
      }
      printf("%-40s %-4.4f K\n", "Info: Input Temperature", sys.T.inKelvin);
    } else if (CheckString(line[0], "Potential")) {
      if (CheckString(line[1], "VDW")) {
        sys.ff.VDW_KIND = sys.ff.VDW_STD_KIND;
        printf("%-40s %-s \n", "Info: Non-truncated potential", "Active");
      } else if (CheckString(line[1], "SHIFT")) {
        sys.ff.VDW_KIND = sys.ff.VDW_SHIFT_KIND;
        printf("%-40s %-s \n", "Info: Shift truncated potential", "Active");
      } else if (CheckString(line[1], "SWITCH")) {
        sys.ff.VDW_KIND = sys.ff.VDW_SWITCH_KIND;
        printf("%-40s %-s \n", "Info: Switch truncated potential", "Active");
      } else if (CheckString(line[1], "EXP6")) {
        sys.ff.VDW_KIND = sys.ff.VDW_EXP6_KIND;
        printf("%-40s %-s \n", "Info: Exp-6 Non-truncated potential", "Active");
      }
    } else if (CheckString(line[0], "LRC")) {
      sys.ff.doTailCorr = checkBool(line[1]);
      if (sys.ff.doTailCorr)
        printf("%-40s %-s \n", "Info: Long Range Correction", "Active");
      else
        printf("%-40s %-s \n", "Info: Long Range Correction", "Inactive");
    } else if (CheckString(line[0], "IPC")) {
      sys.ff.doImpulsePressureCorr = checkBool(line[1]);
      if (sys.ff.doImpulsePressureCorr)
        printf("%-40s %-s \n", "Info: Impulse Pressure Correction", "Active");
      else
        printf("%-40s %-s \n", "Info: Impulse Pressure Correction", "Inactive");
    } else if (CheckString(line[0], "Rswitch")) {
      sys.ff.rswitch = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Switch distance", sys.ff.rswitch);
    } else if (CheckString(line[0], "SubVolumeBox")) {
      if (line.size() == 3) {
        int idx = stringtoi(line[1]);
        uint b = stringtoi(line[2]);
        sys.targetedSwapCollection.AddsubVolumeBox(idx, b);
        sys.intraTargetedSwapCollection.AddsubVolumeBox(idx, b);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 2 values for SubVolumeBox, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeCenter")) {
      if (line.size() >= 5) {
        int idx = stringtoi(line[1]);
        XYZ temp;
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.targetedSwapCollection.AddsubVolumeCenter(idx, temp);
        sys.intraTargetedSwapCollection.AddsubVolumeCenter(idx, temp);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 4 values for SubVolumeCenter, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumePBC")) {
      if (line.size() >= 3) {
        int idx = stringtoi(line[1]);
        sys.targetedSwapCollection.AddsubVolumePBC(idx, line[2]);
        sys.intraTargetedSwapCollection.AddsubVolumePBC(idx, line[2]);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 2 values for SubVolumePBC, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeCenterList")) {
      if (line.size() >= 3) {
        int idx = stringtoi(line[1]);
        std::vector<std::string> temp;
        for (int k = 2; k < (int)line.size(); k++) {
          temp.push_back(line[k]);
        }
        sys.targetedSwapCollection.AddsubVolumeAtomList(idx, temp);
        sys.intraTargetedSwapCollection.AddsubVolumeAtomList(idx, temp);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected atleast 3 values for SubVolumeCenterList, but "
               "received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeDim")) {
      if (line.size() >= 5) {
        int idx = stringtoi(line[1]);
        XYZ temp;
        temp.x = stringtod(line[2]);
        temp.y = stringtod(line[3]);
        temp.z = stringtod(line[4]);
        sys.targetedSwapCollection.AddsubVolumeDimension(idx, temp);
        sys.intraTargetedSwapCollection.AddsubVolumeDimension(idx, temp);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 4 values for SubVolumeDim, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeResidueKind")) {
      if (line.size() >= 3) {
        int idx = stringtoi(line[1]);
        std::vector<std::string> temp;

        for (int k = 2; k < (int)line.size(); k++) {
          temp.push_back(line[k]);
        }
        sys.targetedSwapCollection.AddsubVolumeResKind(idx, temp);
        sys.intraTargetedSwapCollection.AddsubVolumeResKind(idx, temp);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected atleast 2 values for SubVolumeResidueKind, but "
               "received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeRigidSwap")) {
      if (line.size() >= 3) {
        int idx = stringtoi(line[1]);
        bool isRigid = checkBool(line[2]);
        sys.targetedSwapCollection.AddsubVolumeSwapType(idx, isRigid);
        sys.intraTargetedSwapCollection.AddsubVolumeSwapType(idx, isRigid);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 2 values for SubVolumeRigidSwap, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    }
#if ENSEMBLE == GCMC
    else if (CheckString(line[0], "SubVolumeChemPot")) {
      if (line.size() >= 4) {
        int idx = stringtoi(line[1]);
        std::string resName = line[2];
        double value = stringtod(line[3]);
        bool isFugacity = false;
        sys.targetedSwapCollection.AddsubVolumeChemPot(idx, resName, value,
                                                       isFugacity);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 3 values for SubVolumeChemPot, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "SubVolumeFugacity")) {
      if (line.size() >= 4) {
        int idx = stringtoi(line[1]);
        std::string resName = line[2];
        double value = stringtod(line[3]);
        bool isFugacity = true;
        sys.targetedSwapCollection.AddsubVolumeChemPot(idx, resName, value,
                                                       isFugacity);
      } else {
        printf("%-40s %-lu !\n",
               "ERROR: Expected 3 values for SubVolumeFugacity, but received",
               line.size() - 1);
        exit(EXIT_FAILURE);
      }
    }
#endif
    else if (CheckString(line[0], "ExchangeVolumeDim")) {
      if (line.size() == 4) {
        XYZ temp;
        temp.x = stringtod(line[1]);
        temp.y = stringtod(line[2]);
        temp.z = stringtod(line[3]);
        sys.memcVal.subVol = temp;
        sys.intraMemcVal.subVol = temp;
        printf("%-40s %-4.3f %-4.3f %-4.3f A\n",
               "Info: Exchange Sub-Volume Dimensions", temp.x, temp.y, temp.z);
        sys.memcVal.readVol = true;
        sys.intraMemcVal.readVol = true;
      }
    } else if (CheckString(line[0], "ExchangeRatio")) {
      if (line.size() >= 2) {
        printf("%-41s", "Info: ExchangeRatio");
        for (uint i = 1; i < line.size(); i++) {
          uint val = stringtoi(line[i]);
          sys.memcVal.exchangeRatio.push_back(val);
          sys.intraMemcVal.exchangeRatio.push_back(val);
          printf("%-5d", val);
        }
        std::cout << std::endl;
        sys.memcVal.readRatio = true;
        sys.intraMemcVal.readRatio = true;
      }
    } else if (CheckString(line[0], "ExchangeLargeKind")) {
      if (line.size() >= 2) {
        printf("%-41s", "Info: Exchange Large Kind");
        for (uint i = 1; i < line.size(); i++) {
          std::string resName = line[i];
          sys.memcVal.largeKind.push_back(resName);
          sys.intraMemcVal.largeKind.push_back(resName);
          printf("%-5s", resName.c_str());
        }
        std::cout << std::endl;
        sys.memcVal.readLK = true;
        sys.intraMemcVal.readLK = true;
      }
    } else if (CheckString(line[0], "ExchangeSmallKind")) {
      if (line.size() >= 2) {
        printf("%-41s", "Info: Exchange Small Kind");
        for (uint i = 1; i < line.size(); i++) {
          std::string resName = line[i];
          sys.memcVal.smallKind.push_back(resName);
          sys.intraMemcVal.smallKind.push_back(resName);
          printf("%-5s", resName.c_str());
        }
        std::cout << std::endl;
        sys.memcVal.readSK = true;
        sys.intraMemcVal.readSK = true;
      }
    } else if (CheckString(line[0], "SmallKindBackBone")) {
      if ((line.size() % 2) == 0) {
        std::cout << "ERROR: Two atom names must be defined for the backbone "
                     "of each small molecule kind!\n";
        exit(EXIT_FAILURE);
      }
      if (line.size() >= 3) {
        printf("%-41s", "Info: Atom Names in BackBone of Small Molecule Kind ");
        for (uint i = 1; i < line.size() - 1; i += 2) {
          if (i != 1) {
            printf(" , ");
          }
          std::string atom1 = line[i];
          std::string atom2 = line[i + 1];
          sys.memcVal.smallBBAtom1.push_back(atom1);
          sys.memcVal.smallBBAtom2.push_back(atom2);
          sys.intraMemcVal.smallBBAtom1.push_back(atom1);
          sys.intraMemcVal.smallBBAtom2.push_back(atom2);
          printf("%-s - %-s", atom1.c_str(), atom2.c_str());
        }
        std::cout << std::endl;
        sys.memcVal.readSmallBB = true;
        sys.intraMemcVal.readSmallBB = true;
      }
    } else if (CheckString(line[0], "LargeKindBackBone")) {
      if ((line.size() % 2) == 0) {
        std::cout << "ERROR: Two atom names must be defined for the backbone "
                     "of each large molecule kind!\n";
        exit(EXIT_FAILURE);
      }
      if (line.size() >= 3) {
        printf("%-41s", "Info: Atom Names in BackBone of Large Molecule Kind ");
        for (uint i = 1; i < line.size() - 1; i += 2) {
          if (i != 1) {
            printf(" , ");
          }
          std::string atom1 = line[i];
          std::string atom2 = line[i + 1];
          sys.memcVal.largeBBAtom1.push_back(atom1);
          sys.memcVal.largeBBAtom2.push_back(atom2);
          sys.intraMemcVal.largeBBAtom1.push_back(atom1);
          sys.intraMemcVal.largeBBAtom2.push_back(atom2);
          printf("%-s - %-s", atom1.c_str(), atom2.c_str());
        }
        std::cout << std::endl;
        sys.memcVal.readLargeBB = true;
        sys.intraMemcVal.readLargeBB = true;
      }
    } else if (CheckString(line[0], "VDWGeometricSigma")) {
      sys.ff.vdwGeometricSigma = checkBool(line[1]);
      if (sys.ff.vdwGeometricSigma)
        printf("%-40s %-s A\n", "Info: Geometric mean to combine LJ sigma",
               "Active");
    } else if (CheckString(line[0], "Rcut")) {
      sys.ff.cutoff = stringtod(line[1]);
      printf("%-40s %-4.4f A\n", "Info: Cutoff", sys.ff.cutoff);
    } else if (CheckString(line[0], "RcutLow")) {
      sys.ff.cutoffLow = stringtod(line[1]);
      printf("%-40s %-4.4lf A\n", "Info: Short Range Cutoff", sys.ff.cutoffLow);
    } else if (CheckString(line[0], "Exclude")) {
      if (line[1] == sys.exclude.EXC_ONETWO) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETWO_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-TWO");
      } else if (line[1] == sys.exclude.EXC_ONETHREE) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONETHREE_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-THREE");
      } else if (line[1] == sys.exclude.EXC_ONEFOUR) {
        sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
        printf("%-40s %-s \n", "Info: Exclude", "ONE-FOUR");
      }
    } else if (CheckString(line[0], "Ewald")) {
      sys.elect.ewald = checkBool(line[1]);
      sys.elect.readEwald = true;
      if (sys.elect.ewald) {
        printf("%-40s %-s \n", "Info: Ewald Summation", "Active");
      }
    } else if (CheckString(line[0], "ElectroStatic")) {
      sys.elect.enable = checkBool(line[1]);
      sys.elect.readElect = true;
    } else if (CheckString(line[0], "Tolerance")) {
      sys.elect.tolerance = stringtod(line[1]);
      printf("%-40s %-1.3E \n", "Info: Ewald Summation Tolerance",
             sys.elect.tolerance);
    } else if (CheckString(line[0], "RcutCoulomb")) {
      if (line.size() == 3) {
        uint b = stringtoi(line[1]);
        if (b < BOX_TOTAL) {
          sys.elect.cutoffCoulomb[b] = stringtod(line[2]);
          sys.elect.cutoffCoulombRead[b] = true;
          printf("%s %-d %-27s %4.4f A\n", "Info: Box ", b, " CutoffCoulomb",
                 sys.elect.cutoffCoulomb[b]);
        } else {
          std::cout << "ERROR: This simulation requires only " << BOX_TOTAL
                    << " sets of Coulomb Cutoff!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    } else if (CheckString(line[0], "CachedFourier")) {
      sys.elect.cache = checkBool(line[1]);
      sys.elect.readCache = true;
      if (sys.elect.cache) {
        printf("%-40s %-s \n", "Info: Cache Ewald Fourier", "Active");
      } else {
        printf("%-40s %-s \n", "Info: Cache Ewald Fourier", "Inactive");
      }
    } else if (CheckString(line[0], "1-4scaling")) {
      sys.elect.oneFourScale = stringtod(line[1]);
    } else if (CheckString(line[0], "Dielectric")) {
      sys.elect.dielectric = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Dielectric", sys.elect.dielectric);
    } else if (CheckString(line[0], "RunSteps")) {
      sys.step.total = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Total number of steps", sys.step.total);
      if (sys.step.total == 0) {
        in.restart.recalcTrajectory = true;
        printf("%-40s %-s \n", "Info: Recalculate Trajectory", "Active");
        std::cout
            << "ERROR: Recalculate Trajectory is not currently supported!\n";
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "EqSteps")) {
      sys.step.equil = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Number of equilibration steps",
             sys.step.equil);
    } else if (CheckString(line[0], "AdjSteps")) {
      sys.step.adjustment = stringtoi(line[1]);
      printf("%-40s %-lu \n", "Info: Move adjustment frequency",
             sys.step.adjustment);
    } else if (CheckString(line[0], "InitStep")) {
      sys.step.initStep = stringtoi(line[1]);
      sys.step.initStepRead = true;
      printf("%-40s %-lu \n", "Info: InitStep", sys.step.initStep);
    } else if (CheckString(line[0], "PressureCalc")) {
      sys.step.pressureCalc = checkBool(line[1]);
      if (line.size() == 3)
        sys.step.pressureCalcFreq = stringtoi(line[2]);

      if (sys.step.pressureCalc && (line.size() == 2)) {
        std::cout
            << "ERROR: Pressure calculation frequency is not specified!\n";
        exit(EXIT_FAILURE);
      }
      if (!sys.step.pressureCalc)
        printf("%-40s %-s \n", "Info: Pressure calculation", "Inactive");
      else {
        printf("%-40s %-lu \n", "Info: Pressure calculation frequency",
               sys.step.pressureCalcFreq);
      }
    } else if (CheckString(line[0], "ParallelTemperingFreq")) {
      sys.step.parallelTemp = checkBool(line[1]);
      if (line.size() == 3)
        sys.step.parallelTempFreq = stringtoi(line[2]);

      if (sys.step.parallelTemp && (line.size() == 2)) {
        std::cout << "ERROR: Parallel Tempering frequency is not specified!\n";
        exit(EXIT_FAILURE);
      }
      if (!sys.step.parallelTemp)
        printf("%-40s %-s \n", "Info: Parallel Tempering", "Inactive");
      else {
        printf("%-40s %-lu \n", "Info: Parallel Tempering frequency",
               sys.step.parallelTempFreq);
      }
    } else if (CheckString(line[0], "ParallelTemperingAttemptsPerExchange")) {
      sys.step.parallelTemperingAttemptsPerExchange = stringtoi(line[1]);
      printf("%-40s %lu \n", "Info: Number of Attempts Per Exchange Move",
             sys.step.parallelTemperingAttemptsPerExchange);
    } else if (CheckString(line[0], "DisFreq")) {
      sys.moves.displace = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Displacement move frequency",
             sys.moves.displace);
    } else if (CheckString(line[0], "MultiParticleFreq")) {
      sys.moves.multiParticle = stringtod(line[1]);
      if (sys.moves.multiParticle > 0.0) {
        sys.moves.multiParticleEnabled = true;
      }
      printf("%-40s %-4.4f \n", "Info: Multi-Particle move frequency",
             sys.moves.multiParticle);
    } else if (CheckString(line[0], "MultiParticleBrownianFreq")) {
      sys.moves.multiParticleBrownian = stringtod(line[1]);
      if (sys.moves.multiParticleBrownian > 0.0) {
        sys.moves.multiParticleEnabled = true;
      }
      printf("%-40s %-4.4f \n", "Info: Multi-Particle Brownian move frequency",
             sys.moves.multiParticleBrownian);
    } else if (CheckString(line[0], "IntraSwapFreq")) {
      sys.moves.intraSwap = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Intra-Swap move frequency",
             sys.moves.intraSwap);
    } else if (CheckString(line[0], "RegrowthFreq")) {
      sys.moves.regrowth = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Regrowth move frequency",
             sys.moves.regrowth);
    } else if (CheckString(line[0], "CrankShaftFreq")) {
      sys.moves.crankShaft = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Crank-Shaft move frequency",
             sys.moves.crankShaft);
    } else if (CheckString(line[0], "IntraTargetedSwapFreq")) {
      sys.moves.intraTargetedSwap = stringtod(line[1]);
      if (sys.moves.intraTargetedSwap > 0.0) {
        sys.intraTargetedSwapCollection.enable = true;
      }
      printf("%-40s %-4.4f \n", "Info: Intra-Targeted-Swap move frequency",
             sys.moves.intraTargetedSwap);
    } else if (CheckString(line[0], "RotFreq")) {
      sys.moves.rotate = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Rotation move frequency",
             sys.moves.rotate);
    } else if (CheckString(line[0], "IntraMEMC-1Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.intraMemc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: IntraMEMC-2 move frequency",
               sys.moves.intraMemc);
        sys.intraMemcVal.enable = true;
        sys.intraMemcVal.MEMC1 = true;
      }
    } else if (CheckString(line[0], "IntraMEMC-2Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.intraMemc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: IntraMEMC-2 move frequency",
               sys.moves.intraMemc);
        sys.intraMemcVal.enable = true;
        sys.intraMemcVal.MEMC2 = true;
      }
    } else if (CheckString(line[0], "IntraMEMC-3Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.intraMemc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: IntraMEMC-3 move frequency",
               sys.moves.intraMemc);
        sys.intraMemcVal.enable = true;
        sys.intraMemcVal.MEMC3 = true;
      }
    }
#ifdef VARIABLE_VOLUME
    else if (CheckString(line[0], "VolFreq")) {
      sys.moves.volume = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Volume move frequency",
             sys.moves.volume);
    } else if (CheckString(line[0], "useConstantArea")) {
      sys.volume.cstArea = checkBool(line[1]);
      if (sys.volume.cstArea)
        printf("Info: Volume change using constant X-Y area.\n");
      else
        printf("Info: Volume change using constant ratio.\n");
    } else if (CheckString(line[0], "FixVolBox0")) {
      sys.volume.cstVolBox0 = checkBool(line[1]);
      if (sys.volume.cstVolBox0)
        printf("%-40s %-d \n", "Info: Fix volume box", 0);
    }
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
    else if (CheckString(line[0], "SwapFreq")) {
      sys.moves.transfer = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: Molecule swap move frequency",
             sys.moves.transfer);
    } else if (CheckString(line[0], "MEMC-1Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.memc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: MEMC-1 move frequency",
               sys.moves.memc);
        sys.memcVal.enable = true;
        sys.memcVal.MEMC1 = true;
      }
    } else if (CheckString(line[0], "MEMC-2Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.memc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: MEMC-2 move frequency",
               sys.moves.memc);
        sys.memcVal.enable = true;
        sys.memcVal.MEMC2 = true;
      }
    } else if (CheckString(line[0], "MEMC-3Freq")) {
      if (stringtod(line[1]) > 0.0) {
        sys.moves.memc = stringtod(line[1]);
        printf("%-40s %-4.4f \n", "Info: MEMC-3 move frequency",
               sys.moves.memc);
        sys.memcVal.enable = true;
        sys.memcVal.MEMC3 = true;
      }
    } else if (CheckString(line[0], "TargetedSwapFreq")) {
      sys.moves.targetedSwap = stringtod(line[1]);
      if (sys.moves.targetedSwap > 0.0) {
        sys.targetedSwapCollection.enable = true;
      }
      printf("%-40s %-4.4f \n", "Info: Targeted-Swap move frequency",
             sys.moves.targetedSwap);
    } else if (CheckString(line[0], "NeMTMCFreq")) {
      sys.moves.neMolTransfer = stringtod(line[1]);
      printf("%-40s %-4.4f \n", "Info: nonEq Mol-Transfer move frequency",
             sys.moves.neMolTransfer);
      if (sys.moves.neMolTransfer > 0.0) {
        sys.neMTMCVal.enable = true;
      }
    } else if (CheckString(line[0], "LambdaCoulomb")) {
      if (line.size() > 1) {
        sys.neMTMCVal.readLambdaCoulomb = true;
        printf("%-41s", "Info: Lambda Coulomb");
        for (uint i = 1; i < line.size(); i++) {
          double val = stringtod(line[i]);
          sys.neMTMCVal.lambdaCoulomb.push_back(val);
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    } else if (CheckString(line[0], "LambdaVDW")) {
      if (line.size() > 1) {
        sys.neMTMCVal.readLambdaVDW = true;
        printf("%-41s", "Info: Lambda VDW");
        for (uint i = 1; i < line.size(); i++) {
          double val = stringtod(line[i]);
          sys.neMTMCVal.lambdaVDW.push_back(val);
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    } else if (CheckString(line[0], "RelaxingSteps")) {
      if (line.size() > 1) {
        sys.neMTMCVal.readRelaxSteps = true;
        sys.neMTMCVal.relaxSteps = stringtoi(line[1]);
        printf("%-40s %-4d \n", "Info: NeMTMC Relaxing Steps",
               sys.neMTMCVal.relaxSteps);
      }
    } else if (CheckString(line[0], "MultiParticleRelaxing")) {
      if (line.size() > 1) {
        sys.neMTMCVal.MPEnable = checkBool(line[1]);
        sys.neMTMCVal.readMPEnable = true;
        if (sys.neMTMCVal.MPEnable) {
          sys.moves.multiParticleEnabled = sys.neMTMCVal.MPEnable;
          printf("%-40s %s \n", "Info: NeMTMC Relaxing Using MultiParticle",
                 "Active");
        } else {
          printf("%-40s %s \n", "Info: NeMTMC Relaxing Using MultiParticle",
                 "Inactive");
        }
      }
    } else if (CheckString(line[0], "MultiParticleBrownianRelaxing")) {
      if (line.size() > 1) {
        sys.neMTMCVal.MPBEnable = checkBool(line[1]);
        sys.neMTMCVal.readMPBEnable = true;
        if (sys.neMTMCVal.MPBEnable) {
          sys.moves.multiParticleEnabled = sys.neMTMCVal.MPBEnable;
          printf("%-40s %s \n",
                 "Info: NeMTMC Relaxing Using MultiParticleBrownian", "Active");
        } else {
          printf("%-40s %s \n",
                 "Info: NeMTMC Relaxing Using MultiParticleBrownian",
                 "Inactive");
        }
      }
    } else if (CheckString(line[0], "ScalePower")) {
      if (line.size() > 1) {
        sys.neMTMCVal.scalePower = stringtoi(line[1]);
        sys.neMTMCVal.scalePowerRead = true;
        printf("%-40s %-4d \n", "Info: Soft-core scaling power(p)",
               sys.neMTMCVal.scalePower);
      }
    } else if (CheckString(line[0], "ScaleAlpha")) {
      if (line.size() > 1) {
        sys.neMTMCVal.scaleAlpha = stringtod(line[1]);
        sys.neMTMCVal.scaleAlphaRead = true;
        printf("%-40s %-4.4f \n", "Info: Soft-core softness(alpha)",
               sys.neMTMCVal.scaleAlpha);
      }
    } else if (CheckString(line[0], "MinSigma")) {
      if (line.size() > 1) {
        sys.neMTMCVal.scaleSigma = stringtod(line[1]);
        sys.neMTMCVal.scaleSigmaRead = true;
        printf("%-40s %-4.4f A \n", "Info: Soft-core minimum sigma",
               sys.neMTMCVal.scaleSigma);
      }
    } else if (CheckString(line[0], "ScaleCoulomb")) {
      if (line.size() > 1) {
        sys.neMTMCVal.scaleCoulomb = checkBool(line[1]);
        sys.neMTMCVal.scaleCoulombRead = true;
        if (sys.neMTMCVal.scaleCoulomb) {
          printf("%-40s %s \n", "Info: Soft-core for Coulombic interaction",
                 "Active");
        } else {
          printf("%-40s %s \n", "Info: Soft-core for Coulombic interaction",
                 "Inactive");
        }
      }
    } else if (CheckString(line[0], "SampleConfFreq")) {
      if (line.size() > 1) {
        sys.neMTMCVal.conformationProb = stringtod(line[1]);
        sys.neMTMCVal.readConformationProb = true;
        if (sys.neMTMCVal.conformationProb <= 1.0f) {
          printf("%-40s %-4.4f \n",
                 "Info: Intra-Swap/Regrowth Frequency in NeMTMC Relaxing Steps",
                 sys.neMTMCVal.conformationProb);
        }
      }
    } else if (CheckString(line[0], "LambdaVDWLimit")) {
      if (line.size() > 1) {
        sys.neMTMCVal.lambdaLimit = stringtod(line[1]);
        sys.neMTMCVal.readLambdaLimit = true;
        if (sys.neMTMCVal.lambdaLimit <= 1.0f) {
          printf("%-40s %-4.4f \n",
                 "Info: Lambda VDW limit for Intra-Swap move in NeMTMC "
                 "Relaxing Steps",
                 sys.neMTMCVal.lambdaLimit);
        }
      }
    }
#endif
    else if (CheckString(line[0], "CellBasisVector1")) {
      uint box = stringtoi(line[1]);
      if (box < BOX_TOTAL) {
        if (!sys.volume.readCellBasis[box][0]) {
          XYZ temp;
          temp.x = stringtod(line[2]);
          temp.y = stringtod(line[3]);
          temp.z = stringtod(line[4]);
          sys.volume.axis[box].Set(0, temp);
          sys.volume.readCellBasis[box][0] = true;
          sys.volume.hasVolume = sys.volume.ReadCellBasis();
        }
      } else {
        std::cout << "ERROR: This simulation requires only " << BOX_TOTAL
                  << " sets of Cell Basis Vector!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "CellBasisVector2")) {
      uint box = stringtoi(line[1]);
      if (box < BOX_TOTAL) {
        if (!sys.volume.readCellBasis[box][1]) {
          XYZ temp;
          temp.x = stringtod(line[2]);
          temp.y = stringtod(line[3]);
          temp.z = stringtod(line[4]);
          sys.volume.axis[box].Set(1, temp);
          sys.volume.readCellBasis[box][1] = true;
          sys.volume.hasVolume = sys.volume.ReadCellBasis();
        }
      } else {
        std::cout << "ERROR: This simulation requires only " << BOX_TOTAL
                  << " sets of Cell Basis Vector!" << std::endl;
        exit(EXIT_FAILURE);
      }
    } else if (CheckString(line[0], "CellBasisVector3")) {
      uint box = stringtoi(line[1]);
      if (box < BOX_TOTAL) {
        if (!sys.volume.readCellBasis[box][2]) {
          XYZ temp;
          temp.x = stringtod(line[2]);
          temp.y = stringtod(line[3]);
          temp.z = stringtod(line[4]);
          sys.volume.axis[box].Set(2, temp);
          sys.volume.readCellBasis[box][2] = true;
          sys.volume.hasVolume = sys.volume.ReadCellBasis();
        }
      } else {
        std::cout << "ERROR: This simulation requires only " << BOX_TOTAL
                  << " sets of Cell Basis Vector!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if (CheckString(line[0], "CBMC_First")) {
      sys.cbmcTrials.nonbonded.first = stringtoi(line[1]);
      printf("%-40s %-4d \n", "Info: CBMC First atom trials",
             sys.cbmcTrials.nonbonded.first);
    } else if (CheckString(line[0], "CBMC_Nth")) {
      sys.cbmcTrials.nonbonded.nth = stringtoi(line[1]);
      printf("%-40s %-4d \n", "Info: CBMC Secondary atom trials",
             sys.cbmcTrials.nonbonded.nth);
    } else if (CheckString(line[0], "CBMC_Ang")) {
      sys.cbmcTrials.bonded.ang = stringtoi(line[1]);
      printf("%-40s %-4d \n", "Info: CBMC Angle trials",
             sys.cbmcTrials.bonded.ang);
    } else if (CheckString(line[0], "CBMC_Dih")) {
      sys.cbmcTrials.bonded.dih = stringtoi(line[1]);
      printf("%-40s %-4d \n", "Info: CBMC Dihedral trials",
             sys.cbmcTrials.bonded.dih);
    }
#endif
#if ENSEMBLE == GCMC
    else if (CheckString(line[0], "ChemPot")) {
      if (line.size() > 3 && multisim != NULL) {
        std::string resName = line[1];
        if (sys.chemPot.cp.count(resName) > 0) {
          std::cout << "ERROR: Should specify either chemical potential or "
                    << "fugacity at most once for resName " << resName
                    << std::endl;
          exit(EXIT_FAILURE);
        }
        double val = stringtod(line[2 + multisim->worldRank]);
        sys.chemPot.cp[resName] = val;
        printf("%-40s %-6s %-6.4f K\n", "Info: Chemical potential",
               resName.c_str(), val);
      } else if (line.size() != 3) {
        std::cout
            << "ERROR: Chemical potential parameters are not specified!\n";
        exit(EXIT_FAILURE);
      } else {
        std::string resName = line[1];
        if (sys.chemPot.cp.count(resName) > 0) {
          std::cout << "ERROR: Should specify either chemical potential or "
                    << "fugacity at most once for resName " << resName
                    << std::endl;
          exit(EXIT_FAILURE);
        }
        double val = stringtod(line[2]);
        sys.chemPot.cp[resName] = val;
        printf("%-40s %-6s %-6.4f K\n", "Info: Chemical potential",
               resName.c_str(), val);
      }
    } else if (CheckString(line[0], "Fugacity")) {
      if (line.size() != 3) {
        std::cout << "ERROR: Fugacity parameters are not specified!\n";
        exit(EXIT_FAILURE);
      }
      sys.chemPot.isFugacity = true;
      std::string resName = line[1];
      if (sys.chemPot.cp.count(resName) > 0) {
        std::cout << "ERROR: Should specify either chemical potential or "
                  << "fugacity at most once for resName " << resName
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      double val = stringtod(line[2]);
      sys.chemPot.cp[resName] = val * unit::BAR_TO_K_MOLECULE_PER_A3;
      printf("%-40s %-6s %-6.4f bar\n", "Info: Fugacity", resName.c_str(), val);
    }
#endif
#if ENSEMBLE == NVT || ENSEMBLE == NPT
    else if (CheckString(line[0], "LambdaCoulomb")) {
      if (line.size() > 1) {
        sys.freeEn.readLambdaCoulomb = true;
        printf("%-41s", "Info: Lambda Coulomb");
        for (uint i = 1; i < line.size(); i++) {
          double val = stringtod(line[i]);
          sys.freeEn.lambdaCoulomb.push_back(val);
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    } else if (CheckString(line[0], "LambdaVDW")) {
      if (line.size() > 1) {
        sys.freeEn.readLambdaVDW = true;
        printf("%-41s", "Info: Lambda VDW");
        for (uint i = 1; i < line.size(); i++) {
          double val = stringtod(line[i]);
          sys.freeEn.lambdaVDW.push_back(val);
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    } else if (CheckString(line[0], "FreeEnergyCalc")) {
      if (line.size() > 1) {
        sys.freeEn.enable = checkBool(line[1]);
        if (sys.freeEn.enable) {
          printf("%-40s %-s \n", "Info: Free Energy Calculation", "Active");
          if (line.size() > 2) {
            sys.freeEn.frequency = stringtoi(line[2]);
            sys.freeEn.freqRead = true;
            printf("%-40s %-4d \n", "Info: Free Energy Frequency",
                   sys.freeEn.frequency);
          }
        } else {
          printf("%-40s %-s \n", "Info: Free Energy Calculation", "Inactive");
        }
      }
    } else if (CheckString(line[0], "MoleculeType")) {
      if (line.size() > 1) {
        sys.freeEn.molType = line[1];
        sys.freeEn.molTypeRead = true;
        if (line.size() > 2) {
          sys.freeEn.molIndex = stringtoi(line[2]);
          sys.freeEn.molIndexRead = true;
          printf("%-40s %-d in %-s \n", "Info: Free Energy Calc for Molecule",
                 sys.freeEn.molIndex, sys.freeEn.molType.c_str());
        }
      }
    } else if (CheckString(line[0], "InitialState")) {
      if (line.size() > 1) {
        sys.freeEn.iState = stringtoi(line[1]);
        sys.freeEn.iStateRead = true;
        printf("%-40s %-d \n", "Info: Free Energy Calc Lambda state",
               sys.freeEn.iState);
      }
    } else if (CheckString(line[0], "ScalePower")) {
      if (line.size() > 1) {
        sys.freeEn.scalePower = stringtoi(line[1]);
        sys.freeEn.scalePowerRead = true;
        printf("%-40s %-4d \n", "Info: Soft-core scaling power(p)",
               sys.freeEn.scalePower);
      }
    } else if (CheckString(line[0], "ScaleAlpha")) {
      if (line.size() > 1) {
        sys.freeEn.scaleAlpha = stringtod(line[1]);
        sys.freeEn.scaleAlphaRead = true;
        printf("%-40s %-4.4f \n", "Info: Soft-core softness(alpha)",
               sys.freeEn.scaleAlpha);
      }
    } else if (CheckString(line[0], "MinSigma")) {
      if (line.size() > 1) {
        sys.freeEn.scaleSigma = stringtod(line[1]);
        sys.freeEn.scaleSigmaRead = true;
        printf("%-40s %-4.4f \n", "Info: Soft-core minimum sigma",
               sys.freeEn.scaleSigma);
      }
    } else if (CheckString(line[0], "ScaleCoulomb")) {
      if (line.size() > 1) {
        sys.freeEn.scaleCoulomb = checkBool(line[1]);
        sys.freeEn.scaleCoulombRead = true;
        if (sys.freeEn.scaleCoulomb) {
          printf("%-40s %s \n", "Info: Soft-core for Coulombic interaction",
                 "Active");
        } else {
          printf("%-40s %s \n", "Info: Soft-core for Coulombic interaction",
                 "Inactive");
        }
      }
    }
#endif
    else if (CheckString(line[0], "OutputName")) {
      if (multisim != NULL) {
        std::stringstream replicaDirectory;
        replicaDirectory << multisim->replicaOutputDirectoryPath << line[1];
        out.statistics.settings.uniqueStr.val = replicaDirectory.str();
        printf("%-40s %-s \n", "Info: Output name",
               replicaDirectory.str().c_str());
      } else {
        out.statistics.settings.uniqueStr.val = line[1];
        printf("%-40s %-s \n", "Info: Output name", line[1].c_str());
      }
    } else if (CheckString(line[0], "CoordinatesFreq")) {
      out.state.settings.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.state.settings.frequency = stringtoi(line[2]);

      if (out.state.settings.enable && (line.size() == 2))
        out.state.settings.frequency = (ulong)sys.step.total / 10;

      if (out.state.settings.enable) {
        printf("%-40s %-lu \n", "Info: Coordinate frequency",
               out.state.settings.frequency);
      } else
        printf("%-40s %-s \n", "Info: Printing coordinate", "Inactive");
    } else if (CheckString(line[0], "RestartFreq")) {
      out.restart.settings.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.restart.settings.frequency = stringtoi(line[2]);

      if (out.restart.settings.enable && (line.size() == 2))
        out.restart.settings.frequency = (ulong)sys.step.total;

      if (out.restart.settings.enable) {
        printf("%-40s %-lu \n", "Info: Restart frequency",
               out.restart.settings.frequency);
      } else
        printf("%-40s %-s \n", "Info: Printing restart coordinate", "Inactive");
    } else if (CheckString(line[0], "CheckpointFreq")) {
      out.checkpoint.enable = checkBool(line[1]);
      if (line.size() == 3) {
        out.checkpoint.frequency = stringtoi(line[2]);
      }
      if (out.checkpoint.enable) {
        printf("%-40s %-lu \n", "Info: Checkpoint frequency",
               out.checkpoint.frequency);
      } else {
        printf("%-40s %-s \n", "Info: Saving checkpoint", "Inactive");
      }
    } else if (CheckString(line[0], "DCDFreq")) {
      out.state_dcd.settings.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.state_dcd.settings.frequency = stringtoi(line[2]);

      if (out.state_dcd.settings.enable && (line.size() == 2))
        out.state_dcd.settings.frequency = (ulong)sys.step.total / 10;

      if (out.state_dcd.settings.enable) {
        printf("%-40s %-lu \n", "Info: DCD frequency",
               out.state_dcd.settings.frequency);
      } else
        printf("%-40s %-s \n", "Info: Printing DCD ", "Inactive");
    } else if (CheckString(line[0], "ConsoleFreq")) {
      out.console.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.console.frequency = stringtoi(line[2]);

      if (out.console.enable && (line.size() == 2)) {
        if (sys.step.total > 1000) {
          out.console.frequency = (ulong)sys.step.total / 1000;
        } else {
          out.console.frequency = (ulong)sys.step.total / 100;
        }
      }
      if (out.console.enable) {
        printf("%-40s %-lu \n", "Info: Console output frequency",
               out.console.frequency);
      } else
        printf("%-40s %-s \n", "Info: Console output", "Inactive");
    } else if (CheckString(line[0], "BlockAverageFreq")) {
      out.statistics.settings.block.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.statistics.settings.block.frequency = stringtoi(line[2]);

      if (out.statistics.settings.block.enable && (line.size() == 2))
        out.statistics.settings.block.frequency = (ulong)sys.step.total / 100;

      if (out.statistics.settings.block.enable) {
        printf("%-40s %-lu \n", "Info: Average output frequency",
               out.statistics.settings.block.frequency);
      } else
        printf("%-40s %-s \n", "Info: Average output", "Inactive");
    }
#if ENSEMBLE == GCMC
    else if (CheckString(line[0], "HistogramFreq")) {
      out.statistics.settings.hist.enable = checkBool(line[1]);
      if (line.size() == 3)
        out.statistics.settings.hist.frequency = stringtoi(line[2]);

      if (out.statistics.settings.hist.enable && (line.size() == 2)) {
        if (sys.step.total > 1000) {
          out.statistics.settings.hist.frequency = (ulong)sys.step.total / 1000;
        } else {
          out.statistics.settings.hist.frequency = (ulong)sys.step.total / 100;
        }
      }

      if (out.statistics.settings.hist.enable) {
        printf("%-40s %-lu \n", "Info: Histogram output frequency",
               out.statistics.settings.hist.frequency);
      } else
        printf("%-40s %-s \n", "Info: Histogram output", "Inactive");
    } else if (CheckString(line[0], "DistName")) {
      out.state.files.hist.histName = line[1];
    } else if (CheckString(line[0], "HistName")) {
      out.state.files.hist.sampleName = line[1];
    } else if (CheckString(line[0], "RunNumber")) {
      out.state.files.hist.number = line[1];
    } else if (CheckString(line[0], "RunLetter")) {
      out.state.files.hist.letter = line[1];
    } else if (CheckString(line[0], "SampleFreq")) {
      out.state.files.hist.stepsPerHistSample = stringtoi(line[1]);
      printf("%-40s %-d \n", "Info: Histogram sample frequency",
             out.state.files.hist.stepsPerHistSample);
    }
#endif
    else if (CheckString(line[0], "OutEnergy")) {
      out.statistics.vars.energy.block = checkBool(line[1]);
      out.statistics.vars.energy.fluct = checkBool(line[2]);
    } else if (CheckString(line[0], "OutPressure")) {
      out.statistics.vars.pressure.block = checkBool(line[1]);
      out.statistics.vars.pressure.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_PARTICLE_NUMBER
    else if (CheckString(line[0], "OutMolNum")) {
      out.statistics.vars.molNum.block = checkBool(line[1]);
      out.statistics.vars.molNum.fluct = checkBool(line[2]);
    }
#endif
    else if (CheckString(line[0], "OutDensity")) {
      out.statistics.vars.density.block = checkBool(line[1]);
      out.statistics.vars.density.fluct = checkBool(line[2]);
    } else if (CheckString(line[0], "OutSurfaceTension")) {
      out.statistics.vars.surfaceTension.block = checkBool(line[1]);
      out.statistics.vars.surfaceTension.fluct = checkBool(line[2]);
    }
#ifdef VARIABLE_VOLUME
    else if (CheckString(line[0], "OutVolume")) {
      out.statistics.vars.volume.block = checkBool(line[1]);
      out.statistics.vars.volume.fluct = checkBool(line[2]);
    }
#endif
    else if (CheckString(line[0], "Random_Seed")) {
      if (line.size() > 2 && multisim != NULL) {
        in.prng.seed = stringtoi(line[multisim->worldRank + 1]);
      } else {
        in.prng.seed = stringtoi(line[1]);
      }
      if ("INTSEED" == in.prng.kind)
        printf("%-40s %-s \n", "Info: Constant seed", "Active");
      else
        printf("Warning: Constant seed set, but will be ignored.\n");
    } else if (CheckString(line[0], "ParallelTempering_Seed")) {
      in.prngParallelTempering.seed = stringtoi(line[1]);
      if ("INTSEED" == in.prngParallelTempering.kind)
        printf("%-40s %-s \n", "Info: Constant Parallel Tempering seed",
               "Active");
      else
        printf("Warning: Constant Parallel Tempering seed set, but will be "
               "ignored.\n");
    } else if (CheckString(line[0], "ExpertMode")) {
      exptMode = checkBool(line[1]);
      printf("%-40s %-s \n", "Info: Expert Mode",
             exptMode ? "Active" : "Inactive");
    } else {
      std::cout << "Warning: Unknown input " << line[0] << "!" << std::endl;
    }
    // Clear and get ready for the next line
    line.clear();
  }

  //*********** Fill in the default values if not specified ***********//
  fillDefaults();

  //*********** Verify inputs ***********//
  verifyInputs();
  printf("%-40s %-s\n\n", "Finished Reading Input File:", fileName);
}

void ConfigSetup::fillDefaults(void) {
  if (sys.elect.ewald == true) {
    sys.elect.enable = true;
  }

  if (sys.moves.rotate == DBL_MAX) {
    sys.moves.rotate = 0.0;
    printf("%-40s %-4.4f \n", "Default: Rotation move frequency",
           sys.moves.rotate);
  }

  if (sys.moves.intraSwap == DBL_MAX) {
    sys.moves.intraSwap = 0.0;
    printf("%-40s %-4.4f \n", "Default: Intra-Swap move frequency",
           sys.moves.intraSwap);
  }

  if (sys.moves.multiParticle == DBL_MAX) {
    sys.moves.multiParticle = 0.0;
    printf("%-40s %-4.4f \n", "Default: Multi-Particle move frequency",
           sys.moves.multiParticle);
  }

  if (sys.moves.multiParticleBrownian == DBL_MAX) {
    sys.moves.multiParticleBrownian = 0.0;
    printf("%-40s %-4.4f \n", "Default: Multi-Particle Brownian move frequency",
           sys.moves.multiParticleBrownian);
  }

  if (sys.moves.intraMemc == DBL_MAX) {
    sys.moves.intraMemc = 0.0;
    printf("%-40s %-4.4f \n", "Default: Intra-MEMC move frequency",
           sys.moves.intraMemc);
  }

  if (sys.moves.regrowth == DBL_MAX) {
    sys.moves.regrowth = 0.0;
    printf("%-40s %-4.4f \n", "Default: Regrowth move frequency",
           sys.moves.regrowth);
  }

  if (sys.moves.crankShaft == DBL_MAX) {
    sys.moves.crankShaft = 0.0;
    printf("%-40s %-4.4f \n", "Default: Crank-Shaft move frequency",
           sys.moves.crankShaft);
  }

  if (sys.moves.intraTargetedSwap == DBL_MAX) {
    sys.moves.intraTargetedSwap = 0.0;
    printf("%-40s %-4.4f \n", "Default: Targeted-Swap move frequency",
           sys.moves.intraTargetedSwap);
  }

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  if (sys.moves.memc == DBL_MAX) {
    sys.moves.memc = 0.0;
    printf("%-40s %-4.4f \n", "Default: MEMC move frequency", sys.moves.memc);
  }

  if (sys.moves.neMolTransfer == DBL_MAX) {
    sys.moves.neMolTransfer = 0.0;
    printf("%-40s %-4.4f \n", "Default: nonEq Mol-Transfer move frequency",
           sys.moves.neMolTransfer);
  }

  if (sys.moves.targetedSwap == DBL_MAX) {
    sys.moves.targetedSwap = 0.0;
    printf("%-40s %-4.4f \n", "Default: Targeted-Swap move frequency",
           sys.moves.targetedSwap);
  }

  if (sys.neMTMCVal.enable) {
    if (!sys.elect.enable && !sys.neMTMCVal.readLambdaCoulomb) {
      sys.neMTMCVal.lambdaCoulomb.resize(sys.neMTMCVal.lambdaVDW.size(), 0.0);
      sys.neMTMCVal.readLambdaCoulomb = true;
      printf("%-41s", "Default: Lambda Coulomb");
      for (uint i = 0; i < sys.neMTMCVal.lambdaCoulomb.size(); i++) {
        double val = sys.neMTMCVal.lambdaCoulomb[i];
        printf("%-6.3f", val);
      }
      std::cout << std::endl;
    }

    if (sys.neMTMCVal.readLambdaVDW) {
      if (sys.elect.enable && !sys.neMTMCVal.readLambdaCoulomb) {
        sys.neMTMCVal.lambdaCoulomb = sys.neMTMCVal.lambdaVDW;
        sys.neMTMCVal.readLambdaCoulomb = true;
        printf("%-41s", "Default: Lambda Coulomb");
        for (uint i = 0; i < sys.neMTMCVal.lambdaCoulomb.size(); i++) {
          double val = sys.neMTMCVal.lambdaCoulomb[i];
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    }

    if (!sys.neMTMCVal.readMPEnable) {
      sys.neMTMCVal.readMPEnable = true;
      sys.neMTMCVal.MPEnable = false;
      printf("%-40s %s \n", "Default: NeMTMC Relaxing using MultiParticle",
             "Inactive");
    }

    if (!sys.neMTMCVal.readMPBEnable) {
      sys.neMTMCVal.readMPBEnable = true;
      sys.neMTMCVal.MPBEnable = false;
      printf("%-40s %s \n",
             "Default: NeMTMC Relaxing using MultiParticleBrownian",
             "Inactive");
    }

    if (!sys.neMTMCVal.scalePowerRead) {
      sys.neMTMCVal.scalePowerRead = true;
      sys.neMTMCVal.scalePower = 2;
      printf("%-40s %-4d \n", "Default: Soft-core scale power(p)",
             sys.neMTMCVal.scalePower);
    }
    if (!sys.neMTMCVal.scaleAlphaRead) {
      sys.neMTMCVal.scaleAlphaRead = true;
      sys.neMTMCVal.scaleAlpha = 0.5;
      printf("%-40s %-4.4f \n", "Default: Soft-core softness(alpha)",
             sys.neMTMCVal.scaleAlpha);
    }
    if (!sys.neMTMCVal.scaleSigmaRead) {
      sys.neMTMCVal.scaleSigmaRead = true;
      sys.neMTMCVal.scaleSigma = 3.0;
      printf("%-40s %-4.4f A \n", "Default: Soft-core minimum sigma",
             sys.neMTMCVal.scaleSigma);
    }
    if (!sys.neMTMCVal.scaleCoulombRead) {
      sys.neMTMCVal.scaleCoulombRead = true;
      sys.neMTMCVal.scaleCoulomb = false;
      printf("%-40s %s A \n", "Default: Soft-core for Coulombic interaction",
             "Inactive");
    }
    if (!sys.neMTMCVal.readConformationProb) {
      sys.neMTMCVal.readConformationProb = true;
      sys.neMTMCVal.conformationProb = 0.1;
      printf("%-40s %-4.4f A \n",
             "Default: Intra-Swap/Regrowth Frequency in NeMTMC Relaxing Steps",
             sys.neMTMCVal.conformationProb);
    }
    if (!sys.neMTMCVal.readLambdaLimit) {
      sys.neMTMCVal.readLambdaLimit = true;
      sys.neMTMCVal.lambdaLimit = 0.1;
      printf("%-40s %-4.4f A \n",
             "Default: Lambda VDW limit for Intra-Swap move in NeMTMC Relaxing "
             "Steps",
             sys.neMTMCVal.lambdaLimit);
    }
  }

#endif

#if ENSEMBLE == NVT || ENSEMBLE == NPT
  if (sys.freeEn.enable) {
    if (!sys.elect.enable && !sys.freeEn.readLambdaCoulomb) {
      sys.freeEn.lambdaCoulomb.resize(sys.freeEn.lambdaVDW.size(), 0.0);
      sys.freeEn.readLambdaCoulomb = true;
      printf("%-41s", "Default: Lambda Coulomb");
      for (uint i = 0; i < sys.freeEn.lambdaCoulomb.size(); i++) {
        double val = sys.freeEn.lambdaCoulomb[i];
        printf("%-6.3f", val);
      }
      std::cout << std::endl;
    }

    if (sys.freeEn.readLambdaVDW) {
      if (sys.elect.enable && !sys.freeEn.readLambdaCoulomb) {
        sys.freeEn.lambdaCoulomb = sys.freeEn.lambdaVDW;
        sys.freeEn.readLambdaCoulomb = true;
        printf("%-41s", "Default: Lambda Coulomb");
        for (uint i = 0; i < sys.freeEn.lambdaCoulomb.size(); i++) {
          double val = sys.freeEn.lambdaCoulomb[i];
          printf("%-6.3f", val);
        }
        std::cout << std::endl;
      }
    }

    if (!sys.freeEn.scalePowerRead) {
      sys.freeEn.scalePower = 2;
      printf("%-40s %-4d \n", "Default: Soft-core scale power(p)",
             sys.freeEn.scalePower);
    }
    if (!sys.freeEn.scaleAlphaRead) {
      sys.freeEn.scaleAlpha = 0.5;
      printf("%-40s %-4.4f \n", "Default: Soft-core softness(alpha)",
             sys.freeEn.scaleAlpha);
    }
    if (!sys.freeEn.scaleSigmaRead) {
      sys.freeEn.scaleSigma = 3.0;
      printf("%-40s %-4.4f A \n", "Default: Soft-core minimum sigma",
             sys.freeEn.scaleSigma);
    }
    if (!sys.freeEn.scaleCoulombRead) {
      sys.freeEn.scaleCoulomb = false;
      printf("%-40s %s A \n", "Default: Soft-core for Coulombic interaction",
             "Inactive");
    }
  }
#endif

  if (sys.exclude.EXCLUDE_KIND == UINT_MAX) {
    sys.exclude.EXCLUDE_KIND = sys.exclude.EXC_ONEFOUR_KIND;
    printf("%-40s %-s \n", "Default: Exclude", "ONE-FOUR");
  }

  if (sys.elect.oneFourScale == DBL_MAX) {
    if (sys.elect.enable) {
      sys.elect.oneFourScale = 0.0;
      printf("%-40s %-lf \n", "Default: Modified 1-4 Electrostatic scaling",
             sys.elect.oneFourScale);
    }
  }

  if (in.prng.kind == "") {
    in.prng.kind = in.prng.KIND_RANDOM;
    printf("%-40s %-s \n", "Default: Random seed", "Active");
  }

#if ENSEMBLE == GEMC
  if (sys.gemc.kind == UINT_MAX) {
    sys.gemc.kind = mv::GEMC_NVT;
    printf("Default: Running NVT_GEMC\n");
  }
#endif

  if (sys.elect.ewald == true && sys.elect.readCache == false) {
    sys.elect.cache = false;
    printf("%-40s %-s \n", "Default: Cache Ewald Fourier", "Inactive");
  }

  if (sys.elect.ewald == false && sys.elect.cache == true) {
    printf("Warning: Cache Ewald Fourier set, but will be ignored: Ewald "
           "method off.\n");
  }

  if (sys.elect.enable && sys.elect.dielectric == DBL_MAX &&
      in.ffKind.isMARTINI) {
    sys.elect.dielectric = 15.0f;
    printf("%-40s %-4.4f \n", "Default: Dielectric", sys.elect.dielectric);
  }

  if (sys.elect.enable) {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (!sys.elect.cutoffCoulombRead[b]) {
        sys.elect.cutoffCoulomb[b] = sys.ff.cutoff;
        sys.elect.cutoffCoulombRead[b] = true;
        printf("%s %-d %-24s %4.4f A\n", "Default: Box ", b, " CutoffCoulomb",
               sys.elect.cutoffCoulomb[b]);
      }
    }
  }

  if (sys.ff.cutoffLow == DBL_MAX) {
    sys.ff.cutoffLow = 0.0;
    printf("%-40s %-4.4lf \n", "Default: Short Range Cutoff", sys.ff.cutoffLow);
  }

  if (out.statistics.settings.block.enable && in.restart.recalcTrajectory) {
    out.statistics.settings.block.enable = false;
    printf("%-40s \n",
           "Warning: Average output is activated but it will be ignored.");
  }

  if (out.restart.settings.enable && in.restart.recalcTrajectory) {
    out.restart.settings.enable = false;
    printf("%-40s \n", "Warning: Printing restart coordinate is activated but "
                       "it will be ignored.");
  }

  if (out.state.settings.enable && in.restart.recalcTrajectory) {
    out.state.settings.enable = false;
    printf("%-40s \n",
           "Warning: Printing coordinate is activated but it will be ignored.");
  }

  if (out.state_dcd.settings.enable && in.restart.recalcTrajectory) {
    out.state_dcd.settings.enable = false;
    printf("%-40s \n", "Warning: Printing DCD coordinate is activated but it "
                       "will be ignored.");
  }

  out.state.files.psf.name =
      out.statistics.settings.uniqueStr.val + "_merged.psf";
  for (int i = 0; i < BOX_TOTAL; i++) {
    sstrm::Converter toStr;
    std::string numStr = "";
    toStr << i;
    toStr >> numStr;
    out.state.files.pdb.name[i] =
        out.statistics.settings.uniqueStr.val + "_BOX_" + numStr + ".pdb";
    out.restart.files.pdb.name[i] = out.statistics.settings.uniqueStr.val +
                                    "_BOX_" + numStr + "_restart.pdb";
    out.state_dcd.files.dcd.name[i] =
        out.statistics.settings.uniqueStr.val + "_BOX_" + numStr + ".dcd";
    out.restart_dcd.files.dcd.name[i] = out.statistics.settings.uniqueStr.val +
                                        "_BOX_" + numStr + "_restart.coor";
    out.restart_vel.files.dcd.name[i] = out.statistics.settings.uniqueStr.val +
                                        "_BOX_" + numStr + "_restart.vel";
    out.state.files.splitPSF.name[i] =
        out.statistics.settings.uniqueStr.val + "_BOX_" + numStr + ".psf";
  }
  out.state.files.seed.name = out.statistics.settings.uniqueStr.val + ".dat";
}

void ConfigSetup::verifyInputs(void) {
  int i;

#ifdef VARIABLE_PARTICLE_NUMBER
  if (sys.targetedSwapCollection.enable) {
    for (i = 0; i < (int)sys.targetedSwapCollection.targetedSwap.size(); i++) {
      // make sure all required parameter has been set
      sys.targetedSwapCollection.targetedSwap[i].VerifyParm();
    }
  }
#endif

  if (sys.intraTargetedSwapCollection.enable) {
    for (i = 0; i < (int)sys.intraTargetedSwapCollection.targetedSwap.size();
         i++) {
      // make sure all required parameter has been set
      sys.intraTargetedSwapCollection.targetedSwap[i].VerifyParm();
    }
  }

  if (std::abs(sys.moves.multiParticle) > 0.0000001 &&
      std::abs(sys.moves.multiParticleBrownian) > 0.0000001) {
    std::cout << "ERROR: Both multi-Particle and multi-Particle Brownian! "
              << " cannot be used at the same time!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!sys.elect.enable && sys.elect.oneFourScale != DBL_MAX) {
    printf("Warning: 1-4 Electrostatic scaling set, but will be ignored.\n");
    sys.elect.oneFourScale = 0.0;
  }

  if (sys.elect.ewald == false && sys.elect.enable == true) {
    printf("%-40s %-s \n",
           "Warning: Electrostatic calculation with Ewald method", "Inactive");
  }

  if (in.restart.enable && sys.volume.hasVolume) {
    printf(
        "Warning: Cell dimension set, but will be ignored in restart mode.\n");
  }

  if (in.restart.restartFromCheckpoint && !in.restart.enable) {
    std::cout << "ERROR: Restarting from checkpoint file requires"
              << " Restart true!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (in.restart.restartFromBinaryCoorFile && !in.restart.enable) {
    std::cout << "ERROR: Restarting from binary coordinate file(s) requires"
              << " Restart true!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (in.restart.restartFromBinaryVelFile && !in.restart.enable) {
    std::cout << "ERROR: Restarting from binary velocity file(s) requires"
              << " Restart true!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (in.restart.restartFromXSCFile && !in.restart.enable) {
    std::cout << "ERROR: Restarting from extended system file(s) requires"
              << " Restart true!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (in.prng.kind == "RANDOM" && in.prng.seed != UINT_MAX) {
    printf("Warning: Seed value set, but will be ignored.\n");
  }

  // Set output files
  if (out.statistics.settings.uniqueStr.val == "") {
    std::cout << "ERROR: Output name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }

#if ENSEMBLE == GEMC
  if (sys.gemc.kind == mv::GEMC_NPT && sys.gemc.pressure == DBL_MAX) {
    std::cout << "ERROR: Pressure is not specified for NPT-GEMC!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.gemc.kind == mv::GEMC_NVT && sys.gemc.pressure != DBL_MAX) {
    std::cout << "Warning: Input pressure set, but will be ignored in NVT-GEMC"
              << std::endl;
  }

  if (sys.volume.cstVolBox0 && sys.gemc.kind == mv::GEMC_NVT) {
    std::cout << "ERROR: Fix volume of box 0 cannot be applied to NVT-GEMC!\n";
    exit(EXIT_FAILURE);
  }
#endif
#if ENSEMBLE == NPT
  if (sys.gemc.pressure == DBL_MAX) {
    std::cout << "ERROR: Pressure is not specified for NPT!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.volume.cstVolBox0) {
    std::cout
        << "ERROR: Fix volume of box 0 cannot be applied for NPT simulation!\n";
    exit(EXIT_FAILURE);
  }
#endif

#if GOMC_LIB_MPI
  if (sys.step.parallelTemp && in.prngParallelTempering.kind != "INTSEED") {
    std::cout << "ERROR: INTSEED required for parallel tempering!" << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  if (in.prng.kind == "INTSEED" && in.prng.seed == UINT_MAX) {
    std::cout << "ERROR: Seed value is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (in.ffKind.numOfKinds == 0) {
    std::cout << "ERROR: Force field type is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (in.ffKind.numOfKinds > 1) {
    std::cout << "ERROR: Multiple Parameter types are specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if ((!in.ffKind.isMARTINI && !in.ffKind.isEXOTIC) &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND)) {
    std::cout << "Warning: Exclude 1-2 is set for CHARMM type parameter.\n";
  }
  if (in.ffKind.isEXOTIC &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETWO_KIND)) {
    std::cout << "Warning: Exclude 1-2 is set for EXOTIC type parameter.\n";
  }
  if (in.ffKind.isEXOTIC &&
      (sys.exclude.EXCLUDE_KIND == sys.exclude.EXC_ONETHREE_KIND)) {
    std::cout << "Warning: Exclude 1-3 is set for EXOTIC type parameter.\n";
  }
  if (!in.files.param.size()) {
    std::cout << "ERROR: Parameter file name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.ff.VDW_KIND == UINT_MAX) {
    std::cout << "ERROR: Potential type is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (((sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_EXP6_KIND)) &&
      (sys.ff.doTailCorr == false)) {
    std::cout << "Warning: Long Range Correction is Inactive for "
              << "Non-truncated potential." << std::endl;
  }
  if (((sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND)) &&
      sys.ff.doTailCorr) {
    std::cout << "ERROR: LJ long Range Correction is not supported for "
              << "SHIFT or SWITCH potentials." << std::endl;
    exit(EXIT_FAILURE);
  }
  if (((sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_EXP6_KIND)) &&
      (sys.ff.doImpulsePressureCorr && sys.ff.doTailCorr)) {
    std::cout << "ERROR: Impulse Pressure Correction cannot be "
              << "used with LJ long-range corrections." << std::endl;
    exit(EXIT_FAILURE);

  }
  if (((sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND)) &&
      sys.ff.doImpulsePressureCorr) {
    std::cout << "ERROR: Impulse Pressure Correction is not supported "
              << "for SWITCH or SHIFT potentials." << std::endl;
    exit(EXIT_FAILURE);

  }
  if (sys.ff.doImpulsePressureCorr && sys.ff.doTailCorr) {
    std::cout << "ERROR: Both LRC (Long Range Correction) and "
              << "IPC (Impulse Pressure Correction) cannot be used together!"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.ff.cutoff == DBL_MAX) {
    std::cout << "ERROR: Cutoff is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.ff.cutoff < 0.0) {
    std::cout << "ERROR: Potential cutoff cannot be negative!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.ff.cutoffLow < 0.0) {
    sys.ff.cutoffLow = 0.0;
    printf("Warning: RcutLow cannot be negative. Initializing to "
           "zero.\n");
  }
  if (sys.elect.ewald && (sys.elect.tolerance == DBL_MAX)) {
    std::cout << "ERROR: Tolerance is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.step.adjustment == ULONG_MAX) {
    std::cout << "ERROR: Move adjustment frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.step.equil == ULONG_MAX) {
    std::cout << "ERROR: Equilibration steps is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.step.total == ULONG_MAX) {
    std::cout << "ERROR: Total run steps is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.step.adjustment > sys.step.equil && !in.restart.enable &&
      !in.restart.recalcTrajectory) {
    std::cout << "ERROR: Move adjustment frequency cannot be larger "
              << "than the number of equilibration steps!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.step.equil > (sys.step.initStep + sys.step.total) &&
      !in.restart.recalcTrajectory && !in.restart.restartFromCheckpoint) {
    std::cout << "ERROR: The number of equilibration steps cannot exceed "
              << "the total run steps!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (sys.step.equil > (sys.step.initStep + sys.step.total) &&
      !in.restart.recalcTrajectory && !in.restart.restartFromCheckpoint) {
    std::cout << "ERROR: The number of equilibration steps cannot exceed "
              << "the total run steps!" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (sys.moves.displace == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Displacement move frequency is not specified!\n";
      exit(EXIT_FAILURE);
    } else {
      sys.moves.displace = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Displacement move frequency",
             sys.moves.displace);
    }
  }
#if ENSEMBLE == NPT
  if (sys.moves.volume == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Volume move frequency is not specified!"
                << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sys.moves.volume = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Volume move frequency",
             sys.moves.volume);
    }
  }
#endif
#if ENSEMBLE == GEMC
  if (sys.moves.volume == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Volume move frequency is not specified!"
                << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sys.moves.volume = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Volume move frequency",
             sys.moves.volume);
    }
  }

  if (sys.neMTMCVal.enable && sys.moves.transfer == DBL_MAX) {
    sys.moves.transfer = 0.0;
    printf("%-40s %-4.4f \n", "Default: Molecule swap move frequency",
           sys.moves.transfer);
  }

  if (sys.moves.transfer == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Molecule swap move frequency is not specified!"
                << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sys.moves.transfer = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Molecule swap move frequency",
             sys.moves.transfer);
    }
  }

  if (std::abs(sys.moves.displace + sys.moves.rotate + sys.moves.transfer +
               sys.moves.intraSwap + sys.moves.volume + sys.moves.regrowth +
               sys.moves.memc + sys.moves.intraMemc + sys.moves.crankShaft +
               sys.moves.multiParticle + sys.moves.multiParticleBrownian +
               sys.moves.neMolTransfer + sys.moves.targetedSwap +
               sys.moves.intraTargetedSwap - 1.0) > 0.001) {
    std::cout << "ERROR: Sum of move frequencies is not equal to one!\n";
    exit(EXIT_FAILURE);
  }
#elif ENSEMBLE == NPT
  if (sys.moves.volume == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Volume move frequency is not specified!"
                << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sys.moves.volume = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Volume move frequency",
             sys.moves.volume);
    }
  }

  if (std::abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap +
               sys.moves.volume + sys.moves.regrowth + sys.moves.intraMemc +
               sys.moves.crankShaft + sys.moves.multiParticle +
               sys.moves.multiParticleBrownian + sys.moves.intraTargetedSwap -
               1.0) > 0.001) {
    std::cout << "ERROR: Sum of move frequencies is not equal to one!\n";
    exit(EXIT_FAILURE);
  }

#elif ENSEMBLE == GCMC
  if (sys.moves.transfer == DBL_MAX) {
    if (!exptMode) {
      std::cout << "ERROR: Molecule swap move frequency is not specified!"
                << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sys.moves.transfer = 0.0;
      printf("%-40s %-4.4f \n", "ADV USER: Molecule swap move frequency",
             sys.moves.transfer);
    }
  }
  if (std::abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap +
               sys.moves.transfer + sys.moves.regrowth + sys.moves.memc +
               sys.moves.intraMemc + sys.moves.crankShaft +
               sys.moves.multiParticle + sys.moves.multiParticleBrownian +
               sys.moves.neMolTransfer + sys.moves.targetedSwap +
               sys.moves.intraTargetedSwap - 1.0) > 0.001) {
    std::cout << "ERROR: Sum of move frequencies is not equal to one!!\n";
    exit(EXIT_FAILURE);
  }
#else
  if (std::abs(sys.moves.displace + sys.moves.rotate + sys.moves.intraSwap +
               sys.moves.regrowth + sys.moves.intraMemc + sys.moves.crankShaft +
               sys.moves.multiParticle + sys.moves.multiParticleBrownian +
               sys.moves.intraTargetedSwap - 1.0) > 0.001) {
    std::cout << "ERROR: Sum of move frequencies is not equal to one!!\n";
    exit(EXIT_FAILURE);
  }
#endif

  for (i = 0; i < BOX_TOTAL; i++) {
    if (!in.files.pdb.defined[i]) {
      std::cout << "ERROR: PDB file was not specified for box number " << i
                << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  for (i = 0; i < BOX_TOTAL; i++) {
    if (!in.files.psf.defined[i]) {
      std::cout << "ERROR: PSF file was not specified for box number " << i
                << "!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if (in.restart.enable) {
    // ERROR checking to see if if we missed any binary coordinate file
    if (in.restart.restartFromBinaryCoorFile) {
      for (i = 0; i < BOX_TOTAL; i++) {
        if (!in.files.binaryCoorInput.defined[i]) {
          std::cout << "ERROR: Binary coordinate file was not specified for box "
                       "number "
                    << i << "!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    // ERROR checking to see if if we missed any xsc file
    if (in.restart.restartFromXSCFile) {
      for (i = 0; i < BOX_TOTAL; i++) {
        if (!in.files.xscInput.defined[i]) {
          std::cout
              << "ERROR: Extended system file was not specified for box number "
              << i << "!" << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  if (!sys.volume.hasVolume && !in.restart.enable) {
    std::cout << "ERROR: This simulation requires the user define "
              << 3 * BOX_TOTAL << " Cell Basis vectors!" << std::endl;
    for (uint b = 0; b < BOX_TOTAL; b++) {
      for (uint i = 0; i < 3; i++) {
        if (!sys.volume.readCellBasis[b][i]) {
          std::cout << "ERROR: CellBasisVector" << i + 1 << " for Box " << b
                    << " is missing!" << std::endl;
        }
      }
    }
    exit(EXIT_FAILURE);
  }
  if (sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND && sys.ff.rswitch == DBL_MAX) {
    std::cout << "ERROR: Switch distance was not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (((sys.ff.VDW_KIND == sys.ff.VDW_STD_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_SHIFT_KIND) ||
       (sys.ff.VDW_KIND == sys.ff.VDW_EXP6_KIND)) &&
      sys.ff.rswitch != DBL_MAX) {
    std::cout << "Warning: Switch distance set, but will be ignored."
              << std::endl;
  }
  if (sys.ff.VDW_KIND == sys.ff.VDW_SWITCH_KIND &&
      sys.ff.rswitch >= sys.ff.cutoff) {
    std::cout << "ERROR: Switch distance must be less than Rcut!\n";
    exit(EXIT_FAILURE);
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if (sys.cbmcTrials.bonded.ang == UINT_MAX) {
    std::cout << "ERROR: CBMC number of angle trials was not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.cbmcTrials.bonded.dih == UINT_MAX) {
    std::cout << "ERROR: CBMC number of dihedral trials was not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.cbmcTrials.nonbonded.first == UINT_MAX) {
    std::cout << "ERROR: CBMC number of first site trials was not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.cbmcTrials.nonbonded.nth == UINT_MAX) {
    std::cout << "ERROR: CBMC number of nth site trials was not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (sys.memcVal.enable || sys.intraMemcVal.enable) {
    if ((sys.memcVal.MEMC1 && sys.memcVal.MEMC2) ||
        (sys.memcVal.MEMC1 && sys.memcVal.MEMC3) ||
        (sys.memcVal.MEMC2 && sys.memcVal.MEMC3)) {
      std::cout << "ERROR: Multiple MEMC methods were specified, but only one is allowed!\n";
      exit(EXIT_FAILURE);
    }
    if ((sys.intraMemcVal.MEMC1 && sys.intraMemcVal.MEMC2) ||
        (sys.intraMemcVal.MEMC1 && sys.intraMemcVal.MEMC3) ||
        (sys.intraMemcVal.MEMC2 && sys.intraMemcVal.MEMC3)) {
      std::cout << "ERROR: Multiple Intra-MEMC methods are specified, but only one is allowed!\n";
      exit(EXIT_FAILURE);
    }
    if (!sys.memcVal.readVol || !sys.intraMemcVal.readVol) {
      std::cout << "ERROR: In the MEMC method, the Sub-Volume was not specified!\n";
      exit(EXIT_FAILURE);
    }
    if (!sys.memcVal.readRatio || !sys.intraMemcVal.readRatio) {
      std::cout << "ERROR: In the MEMC method, Exchange Ratio was not specified!\n";
      exit(EXIT_FAILURE);
    }
    if (sys.memcVal.largeKind.size() != sys.memcVal.exchangeRatio.size()) {
      std::cout << "ERROR: In the MEMC method, the specified number of Large Kinds was "
                << sys.memcVal.largeKind.size() << ", but "
                << sys.memcVal.exchangeRatio.size()
                << " exchange ratio was specified!\n";
      exit(EXIT_FAILURE);
    }
    if (!sys.memcVal.readSK || !sys.intraMemcVal.readSK) {
      std::cout << "ERROR: In the MEMC method, Small Kind was not specified!\n";
      exit(EXIT_FAILURE);
    }
    if (!sys.memcVal.readLK || !sys.intraMemcVal.readLK) {
      std::cout << "ERROR: In the MEMC method, Large Kind was not specified!\n";
      exit(EXIT_FAILURE);
    }
    if ((sys.memcVal.largeKind.size() != sys.memcVal.smallKind.size()) ||
        (sys.intraMemcVal.largeKind.size() !=
         sys.intraMemcVal.smallKind.size())) {
      std::cout
          << "ERROR: In the MEMC method, the specified number of Large Kinds is not "
          << " equal as specified number of Small Kinds!\n";
      exit(EXIT_FAILURE);
    }
    if (!sys.memcVal.readLargeBB || !sys.intraMemcVal.readLargeBB) {
      std::cout
          << "ERROR: In the MEMC method, Large Kind BackBone was not specified!\n";
      exit(EXIT_FAILURE);
    }
    if (sys.memcVal.largeKind.size() != sys.memcVal.largeBBAtom1.size()) {
      std::cout << "ERROR: In the MEMC method, the specified number of Large Kinds was "
                << sys.memcVal.largeKind.size() << ", but "
                << sys.memcVal.largeBBAtom1.size()
                << " sets of Large Molecule BackBone was specified!\n";
      exit(EXIT_FAILURE);
    }
    if (sys.memcVal.MEMC2 && !sys.memcVal.readSmallBB) {
      std::cout
          << "ERROR: In the MEMC-2 method, Small Kind BackBone was not specified!\n";
      exit(EXIT_FAILURE);
    }

    if (sys.memcVal.MEMC2 &&
        (sys.memcVal.smallKind.size() != sys.memcVal.smallBBAtom1.size())) {
      std::cout
          << "ERROR: In the MEMC-2 method, the specified number of Small Kinds was "
          << sys.memcVal.smallKind.size() << ", but "
          << sys.memcVal.smallBBAtom1.size()
          << " sets of Small Molecule BackBone was specified!\n";
      exit(EXIT_FAILURE);
    }

    if (sys.intraMemcVal.MEMC2 && !sys.intraMemcVal.readSmallBB) {
      std::cout << "ERROR: In the Intra-MEMC-2 method, Small Kind BackBone was not "
                   "specified!\n";
      exit(EXIT_FAILURE);
    }
    if (sys.memcVal.enable && sys.intraMemcVal.enable) {
      if ((sys.memcVal.MEMC1 && !sys.intraMemcVal.MEMC1) ||
          (sys.memcVal.MEMC2 && !sys.intraMemcVal.MEMC2) ||
          (sys.memcVal.MEMC3 && !sys.intraMemcVal.MEMC3)) {
        std::cout << "ERROR: The selected intra-MEMC method was not same as the inter-MEMC method!\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  if (sys.neMTMCVal.enable) {
    if (!sys.neMTMCVal.readLambdaCoulomb) {
      std::cout << "ERROR: Lambda Coulomb states were not defined for "
                << "NeMTMC move! \n";
      exit(EXIT_FAILURE);
    }

    if (!sys.neMTMCVal.readLambdaVDW) {
      std::cout << "ERROR: Lambda VDW states were not defined for "
                << "NeMTMC move! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.neMTMCVal.lambdaCoulomb.size() != sys.neMTMCVal.lambdaVDW.size()) {
      std::cout << "ERROR: Number of Lambda states for VDW and Coulomb "
                << "are not same in NeMTMC move! \n";
      exit(EXIT_FAILURE);
    }

    for (uint i = 0; i < sys.neMTMCVal.lambdaVDW.size(); i++) {
      bool decreasing = false;
      for (uint j = i; j < sys.neMTMCVal.lambdaVDW.size(); j++) {
        if (sys.neMTMCVal.lambdaVDW[i] > sys.neMTMCVal.lambdaVDW[j]) {
          decreasing = true;
        }
      }
      if (decreasing) {
        std::cout << "ERROR: Lambda VDW values are not in increasing order "
                  << "in NeMTMC move! \n";
        exit(EXIT_FAILURE);
      }
    }

    for (uint i = 0; i < sys.neMTMCVal.lambdaCoulomb.size(); i++) {
      bool decreasing = false;
      for (uint j = i; j < sys.neMTMCVal.lambdaCoulomb.size(); j++) {
        if (sys.neMTMCVal.lambdaCoulomb[i] > sys.neMTMCVal.lambdaCoulomb[j]) {
          decreasing = true;
        }
      }
      if (decreasing) {
        std::cout << "ERROR: Lambda Coulomb values are not in increasing "
                  << "order in NeMTMC move! \n";
        exit(EXIT_FAILURE);
      }
    }

    uint last = sys.neMTMCVal.lambdaVDW.size() - 1;
    if (sys.neMTMCVal.lambdaVDW[last] < 0.9999) {
      std::cout << "ERROR: Last Lambda value for VDW is not 1.0 "
                << "in NeMTMC move! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.elect.enable) {
      last = sys.neMTMCVal.lambdaCoulomb.size() - 1;
      if (sys.neMTMCVal.lambdaCoulomb[last] < 0.9999) {
        std::cout << "ERROR: Last Lambda value for Coulomb is not 1.0 "
                  << "in NeMTMC move! \n";
        exit(EXIT_FAILURE);
      }
    }

    if (!sys.neMTMCVal.readRelaxSteps) {
      std::cout << "ERROR: Relaxing steps was not defined for NeMTMC move! \n";
      exit(EXIT_FAILURE);
    } else if (sys.neMTMCVal.relaxSteps == 0) {
      std::cout << "Warning: No thermal relaxing move will be performed in "
                << "NeMTMC move! \n";
    }

    if (sys.neMTMCVal.conformationProb > 1.0f) {
      std::cout
          << "ERROR: Intra-Swap/Regrowth Frequency in NeMTMC Relaxing Steps \n"
          << "       must be less than 1.0! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.neMTMCVal.lambdaLimit > 1.0f) {
      std::cout << "ERROR: Lambda VDW limit for Intra-Swap move in NeMTMC "
                   "Relaxing Steps \n"
                << "       must be less than 1.0! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.moves.multiParticleEnabled) {
      if (sys.neMTMCVal.MPEnable && sys.neMTMCVal.MPBEnable) {
        std::cout
            << "ERROR: Multi-Particle and Multi-Particle Brownian moves cannot "
               "\n"
            << "       be used simultaneously in NeMTMC Relaxing Steps! \n";
        exit(EXIT_FAILURE);
      }
    }

    if (!sys.moves.multiParticleEnabled) {
      if (sys.moves.displace < 1e-7 || sys.moves.rotate < 1e-7) {
        std::cout
            << "ERROR: Displacement and rotation move must be activated in "
            << "NeMTMC Relaxing Steps! \n";
        exit(EXIT_FAILURE);
      }
    }

    if (sys.moves.multiParticleEnabled) {
      if (sys.moves.multiParticle < 1e-7 && sys.neMTMCVal.MPEnable) {
        std::cout << "ERROR: Multi-Particle move must be activated in "
                  << "NeMTMC Relaxing Steps! \n";
        exit(EXIT_FAILURE);
      }
      if (sys.moves.multiParticleBrownian < 1e-7 && sys.neMTMCVal.MPBEnable) {
        std::cout << "ERROR: Multi-Particle Brownian move must be activated in "
                  << "NeMTMC Relaxing Steps! \n";
        exit(EXIT_FAILURE);
      }
    }
  }
#endif
#if ENSEMBLE == NVT || ENSEMBLE == NPT
  if (sys.freeEn.enable) {
    if (sys.ff.cutoffLow > 0.0) {
      sys.ff.cutoffLow = 0.0;
      printf("Warning: Free energy calculations are being used when RcutLow is "
             "not zero (0),\n");
      printf("         which would produce incorrect free energy results.\n");
      printf("         Resetting RcutLow to zero (RcutLow=0) for free energy "
             "calculations!\n");
    }

    if (!sys.freeEn.readLambdaCoulomb) {
      std::cout << "ERROR: Lambda Coulomb states were not defined for "
                << "Free Energy Calculation! \n";
      exit(EXIT_FAILURE);
    }

    if (!sys.freeEn.readLambdaVDW) {
      std::cout << "ERROR: Lambda VDW states were not defined for "
                << "Free Energy Calculation! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.freeEn.lambdaCoulomb.size() != sys.freeEn.lambdaVDW.size()) {
      std::cout << "ERROR: Number of Lambda states for VDW and Coulomb "
                << "are not same in Free Energy Calculation! \n";
      exit(EXIT_FAILURE);
    }

    for (uint i = 0; i < sys.freeEn.lambdaVDW.size(); i++) {
      bool decreasing = false;
      for (uint j = i; j < sys.freeEn.lambdaVDW.size(); j++) {
        if (sys.freeEn.lambdaVDW[i] > sys.freeEn.lambdaVDW[j]) {
          decreasing = true;
        }
      }
      if (decreasing) {
        std::cout << "ERROR: Lambda VDW values are not in increasing order "
                  << "in Free Energy Calculation! \n";
        exit(EXIT_FAILURE);
      }
    }

    for (uint i = 0; i < sys.freeEn.lambdaCoulomb.size(); i++) {
      bool decreasing = false;
      for (uint j = i; j < sys.freeEn.lambdaCoulomb.size(); j++) {
        if (sys.freeEn.lambdaCoulomb[i] > sys.freeEn.lambdaCoulomb[j]) {
          decreasing = true;
        }
      }
      if (decreasing) {
        std::cout << "ERROR: Lambda Coulomb values are not in increasing "
                  << "order in Free Energy Calculation! \n";
        exit(EXIT_FAILURE);
      }
    }

    uint last = sys.freeEn.lambdaVDW.size() - 1;
    if (sys.freeEn.lambdaVDW[last] < 0.9999) {
      std::cout << "ERROR: Last Lambda value for VDW is not 1.0 "
                << "in Free Energy Calculation! \n";
      exit(EXIT_FAILURE);
    }

    if (sys.elect.enable) {
      last = sys.freeEn.lambdaCoulomb.size() - 1;
      if (sys.freeEn.lambdaCoulomb[last] < 0.9999) {
        std::cout << "ERROR: Last Lambda value for Coulomb is not 1.0 "
                  << "in Free Energy Calculation! \n";
        exit(EXIT_FAILURE);
      }
    }

    if (sys.freeEn.lambdaVDW.size() <= sys.freeEn.iState) {
      std::cout << "ERROR: Initial Lambda state is not valid "
                << "in Free Energy Calculation! \n";
      exit(EXIT_FAILURE);
    }

    if (!sys.freeEn.freqRead) {
      std::cout << "ERROR: Frequency of Free Energy Calculation was "
                << "not defined! \n";
      exit(EXIT_FAILURE);
    }
    if (!sys.freeEn.molTypeRead) {
      std::cout << "ERROR: Molecule Type for Free Energy Calculation was "
                << "not defined! \n";
      exit(EXIT_FAILURE);
    }
    if (!sys.freeEn.molIndexRead) {
      std::cout << "ERROR: Molecule Index for Free Energy Calculation was "
                << "not defined! \n";
      exit(EXIT_FAILURE);
    }
#if ENSEMBLE == NVT
    if (sys.step.pressureCalc) {
      if ((sys.freeEn.frequency % sys.step.pressureCalcFreq) != 0) {
        std::cout << "ERROR: Free Energy calculation Freq must be common "
                  << "number of Pressure calculation freq! \n";
        exit(EXIT_FAILURE);
      }
    }
#endif
  }
#endif
  if (sys.T.inKelvin == DBL_MAX) {
    std::cout << "ERROR: Temperature is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.statistics.settings.uniqueStr.val == "") {
    std::cout << "ERROR: Output name is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.console.enable && out.console.frequency == ULONG_MAX) {
    std::cout << "ERROR: Console output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.restart.settings.enable &&
      out.restart.settings.frequency == ULONG_MAX) {
    std::cout << "ERROR: Restart coordinate frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  // If both checkpoint and restart output are enabled,
  // Checkpoint freq must be divisible by restart freq.
  if (out.checkpoint.enable && !(out.restart.settings.enable)) {
    std::cout << "ERROR: CheckpointFreq cannot be used without RestartFreq!\n";
    exit(EXIT_FAILURE);
  }

  if ((out.checkpoint.enable && out.restart.settings.enable) &&
      out.checkpoint.frequency != out.restart.settings.frequency) {
    std::cout << "ERROR: Checkpoint frequency must equal restart frequency!\n";
    std::cout << "Example: RestartFreq 10000; CheckpointFreq 10000\n";
    exit(EXIT_FAILURE);
  }
  if (out.state.settings.enable && out.state.settings.frequency == ULONG_MAX) {
    std::cout << "ERROR: Coordinate frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.state_dcd.settings.enable &&
      out.state_dcd.settings.frequency == ULONG_MAX) {
    std::cout << "ERROR: DCD Coordinate frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }

  if (out.state.settings.enable && out.restart.settings.enable) {
    if (out.state.settings.frequency < out.restart.settings.frequency) {
      if ((out.restart.settings.frequency % out.state.settings.frequency) !=
          0) {
        std::cout
            << "ERROR: Coordinate frequency must be common multiple of \n";
        std::cout << "       restart corrdinate frequency !\n";
        exit(EXIT_FAILURE);
      }
    } else {
      if ((out.state.settings.frequency % out.restart.settings.frequency) !=
          0) {
        std::cout << "ERROR: Restart coordinate frequency must be common "
                     "multiple of \n";
        std::cout << "       corrdinate frequency !\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  if (out.state_dcd.settings.enable && out.restart.settings.enable) {
    if (out.state_dcd.settings.frequency < out.restart.settings.frequency) {
      if ((out.restart.settings.frequency % out.state_dcd.settings.frequency) !=
          0) {
        std::cout << "ERROR: DCD frequency must be common multiple of \n";
        std::cout << "       restart corrdinate frequency !\n";
        exit(EXIT_FAILURE);
      }
    } else {
      if ((out.state_dcd.settings.frequency % out.restart.settings.frequency) !=
          0) {
        std::cout << "ERROR: Restart coordinate frequency must be common "
                     "multiple of \n";
        std::cout << "       DCD frequency !\n";
        exit(EXIT_FAILURE);
      }
    }
  }

  if (out.statistics.settings.block.enable &&
      out.statistics.settings.block.frequency == ULONG_MAX) {
    std::cout << "ERROR: Average output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
#if ENSEMBLE == GCMC
  if (out.statistics.settings.hist.enable &&
      out.statistics.settings.hist.frequency == ULONG_MAX) {
    std::cout << "ERROR: Histogram output frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.state.files.hist.histName == "") {
    std::cout << "ERROR: Distribution file name is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (out.state.files.hist.sampleName == "") {
    std::cout << "ERROR: Histogram file name of is not specified!" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (out.state.files.hist.letter == "") {
    std::cout << "ERROR: Run Letter of histogram file name is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.state.files.hist.number == "") {
    std::cout << "ERROR: Run number of histogram file is not specified!\n";
    exit(EXIT_FAILURE);
  }
  if (out.state.files.hist.stepsPerHistSample == UINT_MAX) {
    std::cout << "ERROR: Histogram output sample frequency is not specified!\n";
    exit(EXIT_FAILURE);
  }
#endif
  if (!out.statistics.settings.block.enable &&
      out.statistics.vars.energy.block) {
    printf("Note: Average output Inactivated. Energy average output will be "
           "ignored.\n");
    out.statistics.vars.energy.block = false;
  }
  if (!out.statistics.settings.block.enable &&
      out.statistics.vars.pressure.block) {
    printf("Note: Average output Inactivated. Pressure average output will be "
           "ignored.\n");
    out.statistics.vars.pressure.block = false;
  }
  if (!sys.step.pressureCalc && out.statistics.vars.pressure.block) {
    printf("Note: Pressure Calculation Inactivated. Pressure average output "
           "will be ignored.\n");
    out.statistics.vars.pressure.block = false;
  }
  if (!sys.step.pressureCalc && out.statistics.vars.surfaceTension.block) {
    printf("Note: Pressure Calculation Inactivated. Surface Tension average "
           "output will be ignored.\n");
    out.statistics.vars.surfaceTension.block = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if (!out.statistics.settings.block.enable &&
      out.statistics.vars.molNum.block) {
    printf("Note: Average output Inactivated. Molecule average output will be "
           "ignored.\n");
    out.statistics.vars.molNum.block = false;
  }
#endif
  if (!out.statistics.settings.block.enable &&
      out.statistics.vars.density.block) {
    printf("Note: Average output Inactivated. Density average output will be "
           "ignored.\n");
    out.statistics.vars.density.block = false;
  }
#ifdef VARIABLE_VOLUME
  if (!out.statistics.settings.block.enable &&
      out.statistics.vars.volume.block) {
    printf("Note: Average output Inactivated. Volume average output will be "
           "ignored.\n");
    out.statistics.vars.volume.block = false;
  }
#endif
  if (!out.console.enable && out.statistics.vars.energy.fluct) {
    printf(
        "Note: Console output Inactivated. Energy output will be ignored.\n");
    out.statistics.vars.energy.fluct = false;
  }
  if (!out.console.enable && out.statistics.vars.pressure.fluct) {
    printf(
        "Note: Console output Inactivated. Pressure output will be ignored.\n");
    out.statistics.vars.pressure.fluct = false;
  }
  if (!sys.step.pressureCalc && out.statistics.vars.pressure.fluct) {
    printf("Note: Pressure Calculation Inactivated. Pressure output will be "
           "ignored.\n");
    out.statistics.vars.pressure.fluct = false;
  }
  if (!sys.step.pressureCalc && out.statistics.vars.surfaceTension.fluct) {
    printf("Note: Pressure Calculation Inactivated. Surface Tension output "
           "will be ignored.\n");
    out.statistics.vars.surfaceTension.fluct = false;
  }
#ifdef VARIABLE_PARTICLE_NUMBER
  if (!out.console.enable && out.statistics.vars.molNum.fluct) {
    printf(
        "Note: Console output Inactivated. Molecule output will be ignored.\n");
  }
#endif
  if (!out.console.enable && out.statistics.vars.density.fluct) {
    printf(
        "Note: Console output Inactivated. Density output will be ignored.\n");
  }
#ifdef VARIABLE_VOLUME
  if (!out.console.enable && out.statistics.vars.volume.fluct) {
    printf(
        "Note: Console output Inactivated. Volume output will be ignored.\n");
  }
#endif
}

const std::string config_setup::PRNGKind::KIND_RANDOM = "RANDOM";
const std::string config_setup::PRNGKind::KIND_SEED = "INTSEED";
const std::string config_setup::FFKind::FF_CHARMM = "CHARMM";
const std::string config_setup::FFKind::FF_EXOTIC = "EXOTIC";
const std::string config_setup::FFKind::FF_MARTINI = "MARTINI";
const std::string config_setup::FFValues::VDW = "VDW";
const std::string config_setup::FFValues::VDW_SHIFT = "VDW_SHIFT";
const std::string config_setup::FFValues::VDW_EXP6 = "VDW_EXP6";
const std::string config_setup::FFValues::VDW_SWITCH = "VDW_SWITCH";
const std::string config_setup::Exclude::EXC_ONETWO = "1-2";
const std::string config_setup::Exclude::EXC_ONETHREE = "1-3";
const std::string config_setup::Exclude::EXC_ONEFOUR = "1-4";

const uint config_setup::FFValues::VDW_STD_KIND = 0;
const uint config_setup::FFValues::VDW_SHIFT_KIND = 1;
const uint config_setup::FFValues::VDW_SWITCH_KIND = 2;
const uint config_setup::FFValues::VDW_EXP6_KIND = 3;
const uint config_setup::Exclude::EXC_ONETWO_KIND = 0;
const uint config_setup::Exclude::EXC_ONETHREE_KIND = 1;
const uint config_setup::Exclude::EXC_ONEFOUR_KIND = 2;
