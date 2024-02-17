/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CONFIG_SETUP_H
#define CONFIG_SETUP_H

#include <iostream> //for cerr, cout;
#include <map>      //for function handle storage.
#include <string>   //for var names, etc.

#include "BasicTypes.h" //for uint, ulong
#include "EnergyTypes.h"
#include "EnsemblePreprocessor.h" //For box total;
#include "FFConst.h"              //For # of param file kinds
#include "InputAbstracts.h"       //For ReadableBase parent class.
#include "InputFileReader.h"      // Input reader
#include "MersenneTwister.h"      //For MTRand::uint32
#include "MoveConst.h"
#include "Reader.h"    //For config reader
#include "UnitConst.h" //For bar --> mol*K / A3 conversion
#include "XYZArray.h"  //For box dimensions.

#if ENSEMBLE == GCMC
#include <sstream> //for reading in variable # of chem. pot.
#endif

#include <sstream> //for prefixing uniqueVal with the pathToReplicaOutputDirectory

#include "GOMC_Config.h" //For PT
#include "ParallelTemperingPreprocessor.h"
#ifdef WIN32
#define OS_SEP '\\'
#else
#define OS_SEP '/'
#endif

namespace config_setup {
/////////////////////////////////////////////////////////////////////////
// Reoccurring structures

// A filename
struct FileName {
  std::string name;
  bool defined;
  FileName(void) {
    defined = false;
    name = "";
  }
  FileName(std::string file, bool def) {
    defined = def;
    name = file;
  }
};

// Multiple filenames
template <uint N> struct FileNames {
  std::string name[N];
  bool defined[N];
  FileNames(void) { std::fill_n(defined, N, false); }
};

/////////////////////////////////////////////////////////////////////////
// Input-specific structures

// Could use "EnableEvent", but restart's enable arguably needs its
// own struct as "frequency" is misleading name for step number
struct RestartSettings {
  bool enable;
  ulong step;
  bool recalcTrajectory;
  bool restartFromCheckpoint;
  bool restartFromBinaryCoorFile;
  bool restartFromBinaryVelFile;
  bool restartFromXSCFile;
  bool operator()(void) { return enable; }
};

// Kinds of Mersenne Twister initialization
struct PRNGKind {
  std::string kind;
  MTRand::uint32 seed;
  bool IsRand(void) const { return str::compare(KIND_RANDOM, kind); }
  bool IsSeed(void) const { return str::compare(KIND_SEED, kind); }
  static const std::string KIND_RANDOM, KIND_SEED;
};

struct FFKind {
  uint numOfKinds;
  bool isCHARMM, isMARTINI, isEXOTIC;
  static const std::string FF_CHARMM, FF_EXOTIC, FF_MARTINI;
};

// Files for input.
struct InFiles {
  std::vector<FileName> param;
  FileNames<BOX_TOTAL> pdb, psf, checkpoint;
  FileNames<BOX_TOTAL> binaryCoorInput, binaryVelInput, xscInput;
  FileName seed;
};

// Input section of config file data.
struct Input {
  RestartSettings restart;
  PRNGKind prng, prngParallelTempering;
  FFKind ffKind;
  InFiles files;
};

/////////////////////////////////////////////////////////////////////////
// System-specific structures

struct Temperature {
  double inKelvin;
};

struct Exclude {
  uint EXCLUDE_KIND;

  static const std::string EXC_ONETWO, EXC_ONETHREE, EXC_ONEFOUR;
  static const uint EXC_ONETWO_KIND, EXC_ONETHREE_KIND, EXC_ONEFOUR_KIND;
};

struct PotentialConfig {
  uint kind;
  double cutoff;
  uint VDW_KIND;
};
struct VDWPot : public PotentialConfig {
  bool doTailCorr;
};
typedef PotentialConfig VDWShift;
struct VDWSwitch : public PotentialConfig {
  double cuton;
};

// Items that effect the system interactions and/or identity, e.g. Temp.
struct FFValues {
  uint VDW_KIND;
  double cutoff, cutoffLow, rswitch;
  bool doTailCorr, vdwGeometricSigma, doImpulsePressureCorr;
  std::string kind;

  static const std::string VDW, VDW_SHIFT, VDW_SWITCH, VDW_EXP6;
  static const uint VDW_STD_KIND, VDW_SHIFT_KIND, VDW_SWITCH_KIND,
      VDW_EXP6_KIND;
};

#if ENSEMBLE == GEMC || ENSEMBLE == NPT

// Items that effect the system interactions and/or identity, e.g. Temp.
struct GEMCKind {
  uint kind;
  double pressure;
};

#endif

struct Step {
  ulong start, total, equil, adjustment, initStep;
  ulong pressureCalcFreq, parallelTempFreq,
      parallelTemperingAttemptsPerExchange;
  bool pressureCalc;
  bool parallelTemp;
  bool initStepRead;
};

// Holds the percentage of each kind of move for this ensemble.
struct MovePercents {
  double displace, rotate, intraSwap, intraMemc, regrowth, crankShaft,
      multiParticle, multiParticleBrownian, intraTargetedSwap;
  bool multiParticleEnabled; // for both multiparticle and multiparticleBrownian
  bool multiParticleLiquid, multiParticleGas; // GEMC: set boxes for MP moves
#ifdef VARIABLE_VOLUME
  double volume;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  double transfer, memc, neMolTransfer, targetedSwap;
#endif
};

struct ElectroStatic {
  bool readEwald;
  bool readWolf;
  bool readDSF;
  bool readElect;
  bool readCache;
  bool enable;
  bool ewald;
  bool wolf;
  bool dsf;
  bool cache;
  bool cutoffCoulombRead[BOX_TOTAL];
  double tolerance;
  double oneFourScale;
  double dielectric;
  double cutoffCoulomb[BOX_TOTAL];
  double wolf_alpha[BOX_TOTAL];
  ElectroStatic(void) {
    std::fill_n(cutoffCoulombRead, BOX_TOTAL, false);
    std::fill_n(cutoffCoulomb, BOX_TOTAL, 0.0);
    std::fill_n(wolf_alpha, BOX_TOTAL, 0.0);
  }
};

struct Volume {
  bool hasVolume, cstArea, cstVolBox0;
  bool readCellBasis[BOX_TOTAL][3];
  XYZArray axis[BOX_TOTAL];
  Volume(void) : hasVolume(false), cstArea(false), cstVolBox0(false) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      axis[b] = XYZArray(3);
      std::fill_n(readCellBasis[b], 3, false);
    }
  }

  bool ReadCellBasis() const {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      for (uint i = 0; i < 3; i++) {
        if (!readCellBasis[b][i])
          return false;
      }
    }
    return true;
  }
};

// If particle number varies (e.g. GCMC, GEMC) load in parameters for
// configurational bias
struct GrowNonbond {
  uint first, nth;
};

// If particle number varies (e.g. GCMC, GEMC) load in parameters for
// configurational bias
struct GrowBond {
  uint ang, dih;
};

struct CBMC {
  GrowNonbond nonbonded;
  GrowBond bonded;
};

#if ENSEMBLE == GCMC
struct ChemicalPotential {
  bool isFugacity;
  std::map<std::string, double> cp;
  ChemicalPotential &operator=(ChemicalPotential const &rhs) {
    isFugacity = rhs.isFugacity;
    cp = rhs.cp;
    return *this;
  }
};
#endif

// Holds information that is required to define the cavity
// for targeted swap
struct TargetSwapParam {
  // defines the center of subVolume
  XYZ subVolumeCenter;
  // defines the dimension of subVolume
  XYZ subVolumeDim;
  // defines the targeted molecule kind
  std::vector<std::string> selectedResKind;
  // defines the atom Indexes to find subVolume center
  std::vector<int> atomList;
  // defines box index for subVolume
  uint selectedBox;
  // defines the subVolume index for error checking
  int subVolumeIdx;
  // defines if we do rigid body insertion/deletion or not
  // defaut value is true
  bool rigid_swap;
  // defines if we have pbc in xyz axis, default is true in all axis
  std::vector<bool> subVolumeBPC;
// specify chemical potential for subVolume
#if ENSEMBLE == GCMC
  ChemicalPotential subVolumeChemPot;
#endif

  bool center_defined, dim_defined, reskind_defined;
  bool box_defined, rigidSwap_defined, atomList_defined;
  bool pbc_defined, chemPot_defined;
  TargetSwapParam(void) {
    center_defined = dim_defined = false;
    reskind_defined = box_defined = false;
    rigidSwap_defined = pbc_defined = false;
    rigid_swap = atomList_defined = false;
    chemPot_defined = false;
    subVolumeIdx = 0;
    subVolumeBPC.resize(3, true);
  }

  TargetSwapParam &operator=(TargetSwapParam const &rhs) {
    subVolumeCenter = rhs.subVolumeCenter;
    subVolumeDim = rhs.subVolumeDim;
    selectedResKind = rhs.selectedResKind;
    atomList = rhs.atomList;
    selectedBox = rhs.selectedBox;
    subVolumeIdx = rhs.subVolumeIdx;
    rigid_swap = rhs.rigid_swap;
    subVolumeBPC = rhs.subVolumeBPC;
#if ENSEMBLE == GCMC
    subVolumeChemPot = rhs.subVolumeChemPot;
#endif
    // copy boolean parameters
    center_defined = rhs.center_defined;
    dim_defined = rhs.dim_defined;
    reskind_defined = rhs.reskind_defined;
    box_defined = rhs.box_defined;
    rigidSwap_defined = rhs.rigidSwap_defined;
    atomList_defined = rhs.atomList_defined;
    pbc_defined = rhs.pbc_defined;
    return *this;
  }

  // Make sure all parameters have been set
  void VerifyParm() {
    bool allSet = true;
    if (!box_defined) {
      printf("Error: Box has not been defined for subVolume index %d!\n",
             subVolumeIdx);
      allSet = false;
    }
    if (!center_defined && !atomList_defined) {
      printf("Error: Center or atom list (to find the center) has not been "
             "defined for subVolume index %d!\n",
             subVolumeIdx);
      allSet = false;
    }
    if (!dim_defined) {
      printf("Error: dimension has not been defined for subVolume index %d!\n",
             subVolumeIdx);
      allSet = false;
    }

    if (subVolumeDim.x < 0.0 || subVolumeDim.y < 0.0 || subVolumeDim.z < 0.0) {
      printf("Error: SubVolume dimension (%g %g %g) cannot be negative for "
             "subVolume index %d!\n",
             subVolumeDim.x, subVolumeDim.y, subVolumeDim.z, subVolumeIdx);
      allSet = false;
    }

    if (!reskind_defined) {
      printf(
          "Error: residue kind has not been defined for subVolume index %d!\n",
          subVolumeIdx);
      allSet = false;
    }
    if (!rigidSwap_defined) {
      printf("Default: Rigid body swap type has been defined for subVolume "
             "index %d!\n",
             subVolumeIdx);
      rigid_swap = true;
    }
    if (!pbc_defined) {
      printf("Default: XYZ PBC has been defined for subVolume index %d!\n",
             subVolumeIdx);
      pbc_defined = true;
    }

#if ENSEMBLE == GCMC
    // make sure the resname exist in the list of targeted reskind
    if (chemPot_defined) {
      bool notExist = false;
      std::map<std::string, double>::const_iterator start =
          subVolumeChemPot.cp.begin();
      std::map<std::string, double>::const_iterator end =
          subVolumeChemPot.cp.end();
      if ((selectedResKind[0] != "ALL") && reskind_defined) {
        while (start != end) {
          std::string resname = start->first;
          notExist = (std::find(selectedResKind.begin(), selectedResKind.end(),
                                resname) == selectedResKind.end());

          if (notExist) {
            allSet = false;
            printf("Error: residue name %s has not been defined in "
                   "SubVolumeResidueKind list for subVolume index %d!\n",
                   resname.c_str(), subVolumeIdx);
          }

          ++start;
        }
      }
    }
#endif

    if (!allSet)
      exit(EXIT_FAILURE);
  }
};

struct TargetSwapCollection {
  TargetSwapCollection(void) { enable = false; }
  // Search for cavIdx in the targetSwap. If it exists,
  // returns true + the index to targetSwap vector
  bool SearchExisting(const int &subVIdx, int &idx) {
    for (int i = 0; i < (int)targetedSwap.size(); i++) {
      if (targetedSwap[i].subVolumeIdx == subVIdx) {
        idx = i;
        return true;
      }
    }
    return false;
  }

  // add a new subVolume
  void AddsubVolumeBox(const int &subVIdx, uint &box) {
    int idx = 0;
    if (box < BOXES_WITH_U_NB) {
      if (!SearchExisting(subVIdx, idx)) {
        // If the subVolume index did not exist, add one
        TargetSwapParam tempPar;
        tempPar.subVolumeIdx = subVIdx;
        tempPar.selectedBox = box;
        tempPar.box_defined = true;
        targetedSwap.push_back(tempPar);
      } else {
        // If subVolume index exist and subvolume box is defined
        if (targetedSwap[idx].box_defined) {
          printf("Error: The subVolume index %d has already been defined for "
                 "Box %d!\n",
                 subVIdx, targetedSwap[idx].selectedBox);
          printf("       Please use different subVolume index.\n");
          exit(EXIT_FAILURE);
        } else {
          targetedSwap[idx].selectedBox = box;
          targetedSwap[idx].box_defined = true;
        }
      }
    } else {
      printf("Error: Subvolume index %d cannot be set for box %d!\n", subVIdx,
             box);
#if ENSEMBLE == GCMC
      printf("       Maximum box index for this simulation is %d!\n",
             BOXES_WITH_U_NB - 1);
#else
      printf("       Maximum box index for this simulation is %d!\n",
             BOXES_WITH_U_NB - 1);
#endif
      exit(EXIT_FAILURE);
    }
  }

  // add dimension of subVolume to subVIdx
  void AddsubVolumeDimension(const int &subVIdx, XYZ &dimension) {
    int idx = 0;
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.subVolumeDim = dimension;
      tempPar.dim_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist and subvolume dimension is defined
      if (targetedSwap[idx].dim_defined) {
        printf("Error: The subVolume dimension (%g, %g, %g) has already been "
               "defined for subVolume index %d!\n",
               targetedSwap[idx].subVolumeDim.x,
               targetedSwap[idx].subVolumeDim.y,
               targetedSwap[idx].subVolumeDim.z, subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else {
        targetedSwap[idx].subVolumeDim = dimension;
        targetedSwap[idx].dim_defined = true;
      }
    }
  }

  // add center subVolume to subVIdx
  void AddsubVolumeCenter(const int &subVIdx, XYZ &center) {
    int idx = 0;
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.subVolumeCenter = center;
      tempPar.center_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist and and subvolume center is defined
      if (targetedSwap[idx].center_defined) {
        printf("Error: The subVolume center (%g, %g, %g) has already been "
               "defined for subVolume index %d!\n",
               targetedSwap[idx].subVolumeCenter.x,
               targetedSwap[idx].subVolumeCenter.y,
               targetedSwap[idx].subVolumeCenter.z, subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else if (targetedSwap[idx].atomList_defined) {
        printf("Error: Atom list has already been defined to find center in "
               "subVolume index %d!\n",
               subVIdx);
        printf("       One of the SubVolumeCenter or SubVolumeCenterList must "
               "be defined.\n");
        exit(EXIT_FAILURE);
      } else {
        targetedSwap[idx].subVolumeCenter = center;
        targetedSwap[idx].center_defined = true;
      }
    }
  }

  // add center subVolume to subVIdx
  void AddsubVolumeAtomList(const int &subVIdx,
                            std::vector<std::string> &atoms) {
    int idx = 0;
    std::vector<int> alist = ParsAtomList(subVIdx, atoms);
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.atomList = alist;
      tempPar.atomList_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist and and subvolume center is defined
      if (targetedSwap[idx].atomList_defined) {
        printf("Error: The atom list has already been defined to find center "
               "for subVolume index %d!\n",
               subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else if (targetedSwap[idx].center_defined) {
        printf("Error: Subvolume center (%g, %g, %g) has been defined for in "
               "subVolume index %d!\n",
               targetedSwap[idx].subVolumeCenter.x,
               targetedSwap[idx].subVolumeCenter.y,
               targetedSwap[idx].subVolumeCenter.z, subVIdx);
        printf("       One of the SubVolumeCenter or SubVolumeCenterList must "
               "be defined.\n");
        exit(EXIT_FAILURE);
      } else {
        targetedSwap[idx].atomList = alist;
        targetedSwap[idx].atomList_defined = true;
      }
    }
  }

  int stringtoui(const std::string &s) {
    int i = std::stoi(s);
    std::string newstr = std::to_string(i);

    if (i < 0) {
      printf("Error: Expected to receive unsigned integer, but received %s!\n",
             s.c_str());
      exit(EXIT_FAILURE);
    }
    if (s != newstr) {
      printf("Warning: Converting %s to %s in parsing integer!\n", s.c_str(),
             newstr.c_str());
    }

    return i;
  }

  // Parse the list of atoms. If user uses a range of atom (e.g. 5-10)
  // This function would handle it properly
  std::vector<int> ParsAtomList(const int &subVIdx,
                                std::vector<std::string> &atoms) {
    const char dash = '-';
    size_t pos = 0;
    int start = 0, end = 0;
    std::vector<int> list;
    for (int i = 0; i < (int)atoms.size(); ++i) {
      if (atoms[i] == "-") {
        printf("Error: In subVolume Index %d, if you are trying to use a range "
               "of atom indices, use '-' without any space!\n",
               subVIdx);
        exit(EXIT_FAILURE);
      }
      // check if string has dash in it
      pos = atoms[i].find(dash);
      if (pos != std::string::npos && pos != 0) {
        start = stringtoui(atoms[i].substr(0, pos));
        // dash has length of 1
        end = stringtoui(atoms[i].substr(pos + 1, atoms[i].length()));
        while (start <= end) {
          list.push_back(start);
          ++start;
        }

      } else {
        // it's single atom index
        list.push_back(stringtoui(atoms[i]));
      }
    }
    return list;
  }

  // add targeted residue kind of subVolume to subVIdx
  void AddsubVolumeResKind(const int &subVIdx, std::vector<std::string> &kind) {
    CheckSelectedKind(subVIdx, kind);
    int idx = 0;
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      // if there is any kind defined!
      if (kind.size()) {
        tempPar.selectedResKind = kind;
        tempPar.reskind_defined = true;
        targetedSwap.push_back(tempPar);
      }
    } else {
      // If subVolume index exist and and reskind is defined
      if (targetedSwap[idx].reskind_defined) {
        printf("Error: The targeted residue kind has already been defined for "
               "subVolume index %d!\n",
               subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else {
        // if there is any kind defined!
        if (kind.size()) {
          targetedSwap[idx].selectedResKind = kind;
          targetedSwap[idx].reskind_defined = true;
        }
      }
    }
  }
  // Check if user defined multiple reskind
  void CheckSelectedKind(const int &subVIdx, std::vector<std::string> &kind) {
    std::vector<std::string> newKind = kind;
    std::vector<std::string>::iterator ip;
    std::sort(newKind.begin(), newKind.end());
    ip = std::unique(newKind.begin(), newKind.end());
    // Find Uniquie reskind
    newKind.resize(std::distance(newKind.begin(), ip));
    if (newKind.size() != kind.size()) {
      printf("Warning: Duplicated residue kind was defined for subVolume index "
             "%d!\n",
             subVIdx);
      printf(
          "Warning: Proceed with unique residue kind for subVolume index %d!\n",
          subVIdx);
    }

    bool selectedAll =
        (std::find(newKind.begin(), newKind.end(), "ALL") != newKind.end());
    selectedAll |=
        (std::find(newKind.begin(), newKind.end(), "all") != newKind.end());
    selectedAll |=
        (std::find(newKind.begin(), newKind.end(), "All") != newKind.end());
    if (selectedAll) {
      if (newKind.size() > 1) {
        printf("Warning: %lu additional residue kinds were defined for "
               "subVolume index %d, while using all residues!\n",
               newKind.size() - 1, subVIdx);
        printf(
            "Warning: Proceed with all residue kinds for subVolume index %d!\n",
            subVIdx);
      }
      newKind.clear();
      newKind.push_back("ALL");
    }

    kind = newKind;
  }

  // add a swapType for subvolume
  void AddsubVolumeSwapType(const int &subVIdx, bool &isRigid) {
    int idx = 0;
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.rigid_swap = isRigid;
      tempPar.rigidSwap_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist and subvolume box is defined
      if (targetedSwap[idx].rigidSwap_defined) {
        printf("Error: The swap type has already been defined for subVolume "
               "index %d!\n",
               subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else {
        targetedSwap[idx].rigid_swap = isRigid;
        targetedSwap[idx].rigidSwap_defined = true;
      }
    }
  }

  // add a PBC mode for subvolume
  void AddsubVolumePBC(const int &subVIdx, const std::string &pbc) {
    int idx = 0;
    std::vector<bool> pbcMode(3, false);
    char upper;
    for (uint k = 0; k < pbc.length(); k++) {
      upper = toupper(pbc[k]);
      if (upper == 'X') {
        pbcMode[0] = true;
      } else if (upper == 'Y') {
        pbcMode[1] = true;
      } else if (upper == 'Z') {
        pbcMode[2] = true;
      } else {
        printf("Error: Unknown option '%c' in '%s' for PBC type in subVolume "
               "index %d!\n",
               pbc[k], pbc.c_str(), subVIdx);
        exit(EXIT_FAILURE);
      }
    }

    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.subVolumeBPC = pbcMode;
      tempPar.pbc_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist and subvolume box is defined
      if (targetedSwap[idx].pbc_defined) {
        printf("Error: The PBC mode has already been defined for subVolume "
               "index %d!\n",
               subVIdx);
        printf("       Please use different subVolume index.\n");
        exit(EXIT_FAILURE);
      } else {
        targetedSwap[idx].subVolumeBPC = pbcMode;
        targetedSwap[idx].pbc_defined = true;
      }
    }
  }

#if ENSEMBLE == GCMC
  // set chemical potential for subvolume
  void AddsubVolumeChemPot(const int &subVIdx, const std::string &resName,
                           const double &cpValue, const bool &isFugacity) {
    int idx = 0;
    double value =
        cpValue * (isFugacity ? unit::BAR_TO_K_MOLECULE_PER_A3 : 1.0);
    if (!SearchExisting(subVIdx, idx)) {
      // If the subVolume index did not exist, add one
      TargetSwapParam tempPar;
      tempPar.subVolumeIdx = subVIdx;
      tempPar.subVolumeChemPot.cp[resName] = value;
      tempPar.subVolumeChemPot.isFugacity = isFugacity;
      tempPar.chemPot_defined = true;
      targetedSwap.push_back(tempPar);
    } else {
      // If subVolume index exist
      targetedSwap[idx].subVolumeChemPot.cp[resName] = value;
      targetedSwap[idx].subVolumeChemPot.isFugacity = isFugacity;
      targetedSwap[idx].chemPot_defined = true;
    }
  }
#endif

public:
  std::vector<TargetSwapParam> targetedSwap;
  bool enable;
};

struct MEMCVal {
  bool enable, readVol, readRatio, readSmallBB, readLargeBB;
  bool readSK, readLK;
  bool MEMC1, MEMC2, MEMC3;
  XYZ subVol;
  std::vector<std::string> smallKind, largeKind;
  std::vector<uint> exchangeRatio;
  std::vector<std::string> smallBBAtom1, smallBBAtom2;
  std::vector<std::string> largeBBAtom1, largeBBAtom2;
  MEMCVal(void) {
    MEMC1 = MEMC2 = MEMC3 = false;
    readVol = readRatio = readSmallBB = false;
    readLargeBB = readSK = readLK = false;
  }
};

struct NEMTMCVal {
  // relaxing parameter
  uint relaxSteps;
  double lambdaLimit, conformationProb;
  bool MPEnable, MPBEnable;
  bool readRelaxSteps, readMPEnable, readMPBEnable;
  bool readLambdaLimit, readConformationProb;
  // scaling parameter
  uint scalePower;
  double scaleAlpha, scaleSigma;
  bool scaleCoulomb;
  bool enable, readLambdaCoulomb, readLambdaVDW;
  bool scalePowerRead, scaleAlphaRead, scaleSigmaRead, scaleCoulombRead;
  std::vector<double> lambdaCoulomb, lambdaVDW;
  NEMTMCVal(void) {
    readLambdaCoulomb = readRelaxSteps = false;
    readMPEnable = MPEnable = readLambdaVDW = enable = false;
    scalePowerRead = scaleAlphaRead = scaleSigmaRead = scaleCoulombRead = false;
    MPBEnable = readMPBEnable = readLambdaLimit = readConformationProb = false;
  }
};

struct FreeEnergy {
  bool enable, readLambdaCoulomb, readLambdaVDW, freqRead;
  bool molTypeRead, molIndexRead, iStateRead;
  uint frequency, molIndex, iState;
  // scaling parameter
  uint scalePower;
  double scaleAlpha, scaleSigma;
  bool scaleCoulomb;
  bool scalePowerRead, scaleAlphaRead, scaleSigmaRead, scaleCoulombRead;
  std::string molType;
  std::vector<double> lambdaCoulomb, lambdaVDW;
  FreeEnergy(void) {
    readLambdaCoulomb = readLambdaVDW = enable = freqRead = false;
    molTypeRead = molIndexRead = iStateRead = false;
    scalePowerRead = scaleAlphaRead = scaleSigmaRead = scaleCoulombRead = false;
  }
};

struct SystemVals {
  ElectroStatic elect;
  Temperature T;
  FFValues ff;
  Exclude exclude;
  Step step;
  MovePercents moves;
  Volume volume; // May go unused
  CBMC cbmcTrials;
  MEMCVal memcVal, intraMemcVal;
  NEMTMCVal neMTMCVal;
  FreeEnergy freeEn;
  TargetSwapCollection targetedSwapCollection;
  TargetSwapCollection intraTargetedSwapCollection;
#if ENSEMBLE == GCMC
  ChemicalPotential chemPot;
#elif ENSEMBLE == GEMC || ENSEMBLE == NPT
  GEMCKind gemc;
#endif
};

/////////////////////////////////////////////////////////////////////////
// Output-specific structures

struct EventSettings { /* : ReadableStepDependentBase*/
  bool enable;
  ulong frequency;
  bool operator()(void) { return enable; }
};

struct UniqueStr { /* : ReadableBase*/
  std::string val;
};

struct HistFiles { /* : ReadableBase*/
  std::string histName, number, letter, sampleName;
  uint stepsPerHistSample;
};

// Files for output.
struct OutFiles {
  /* For split pdb, psf, and dcd files , BOX 0 and BOX 1 */
  FileNames<BOX_TOTAL> pdb, splitPSF, dcd;

  /* For merged PSF */
  FileName psf, seed;
  HistFiles hist;
};
struct Settings {
  EventSettings block, hist;
  UniqueStr uniqueStr;
};

// Enables for each variable that can be tracked
struct OutputEnables {
  bool block, fluct;
};

struct TrackedVars {
  OutputEnables energy, pressure;
#ifdef VARIABLE_VOLUME
  OutputEnables volume;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  OutputEnables molNum;
#endif
  OutputEnables density;
  OutputEnables surfaceTension;
};

struct SysState {
  EventSettings settings;
  OutFiles files;
};
struct Statistics {
  Settings settings;
  TrackedVars vars;
};
struct Output {
  SysState state, restart, state_dcd, restart_dcd;
  SysState restart_vel;
  Statistics statistics;
  EventSettings console, checkpoint;
};

} // namespace config_setup

class ConfigSetup {
public:
  config_setup::Input in;
  config_setup::Output out;
  config_setup::SystemVals sys;
  bool exptMode;
  ConfigSetup(void);
  void Init(const char *fileName, MultiSim const *const &multisim);

private:
  void fillDefaults(void);
  int stringtoi(const std::string &s);
  double stringtod(const std::string &s);
  bool checkBool(std::string str);
  bool isBool(std::string str);
  bool CheckString(std::string str1, std::string str2);
  void verifyInputs(void);
  InputFileReader reader;
};

#endif /*CONFIG_SETUP_H*/
