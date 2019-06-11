/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CONFIGSETUP_H
#define CONFIGSETUP_H

#include <map> //for function handle storage.
#include <iostream> //for cerr, cout;
#include <string> //for var names, etc.

#include "InputAbstracts.h" //For ReadableBase parent class.
#include "InputFileReader.h" // Input reader
#include "BasicTypes.h" //for uint, ulong
#include "EnsemblePreprocessor.h" //For box total;
#include "Reader.h" //For config reader
#include "FFConst.h" //For # of param file kinds
#include "MersenneTwister.h" //For MTRand::uint32
#include "XYZArray.h" //For box dimensions.
#include "MoveConst.h"
#include "UnitConst.h" //For bar --> mol*K / A3 conversion

#if ENSEMBLE == GCMC
#include <sstream>  //for reading in variable # of chem. pot.
#endif

namespace config_setup
{
/////////////////////////////////////////////////////////////////////////
// Reoccurring structures

//A filename
struct FileName {
  std::string name;
};

//Multiple filenames
template <uint N>
struct FileNames {
  std::string name[N];
};

/////////////////////////////////////////////////////////////////////////
// Input-specific structures

//Could use "EnableEvent", but restart's enable arguably needs its
//own struct as "frequency" is misleading name for step number
struct RestartSettings {
  bool enable;
  ulong step;
  bool recalcTrajectory;
  bool restartFromCheckpoint;
  bool operator()(void)
  {
    return enable;
  }
};

//Kinds of Mersenne Twister initialization
struct PRNGKind {
  std::string kind;
  MTRand::uint32 seed;
  bool IsRand(void) const
  {
    return str::compare(KIND_RANDOM, kind);
  }
  bool IsSeed(void) const
  {
    return str::compare(KIND_SEED, kind);
  }
  bool IsRestart(void) const
  {
    return str::compare(KIND_RESTART, kind);
  }
  static const std::string KIND_RANDOM, KIND_SEED, KIND_RESTART;
};

struct FFKind {
  uint numOfKinds;
  bool isCHARMM, isMARTINI, isEXOTIC;
  static const std::string FF_CHARMM, FF_EXOTIC, FF_MARTINI;
};

//Files for input.
struct InFiles {
  FileName param;
  FileNames<BOX_TOTAL> pdb, psf;
  FileName seed;
};

//Input section of config file data.
struct Input {
  RestartSettings restart;
  PRNGKind prng;
  FFKind ffKind;
  InFiles files;
};



/////////////////////////////////////////////////////////////////////////
// System-specific structures

struct Temperature {
  real inKelvin;
};

struct Exclude {
  uint EXCLUDE_KIND;

  static const std::string EXC_ONETWO, EXC_ONETHREE, EXC_ONEFOUR;
  static const uint EXC_ONETWO_KIND, EXC_ONETHREE_KIND, EXC_ONEFOUR_KIND;
};

struct PotentialConfig {
  uint kind;
  real cutoff;
  uint VDW_KIND;
};
struct VDWPot : public PotentialConfig {
  bool doTailCorr;
};
typedef PotentialConfig VDWShift;
struct VDWSwitch : public PotentialConfig {
  real cuton;
};

//Items that effect the system interactions and/or identity, e.g. Temp.
struct FFValues {
  uint VDW_KIND;
  real cutoff, cutoffLow, rswitch;
  bool doTailCorr, vdwGeometricSigma;
  std::string kind;

  static const std::string VDW, VDW_SHIFT, VDW_SWITCH;
  static const uint VDW_STD_KIND, VDW_SHIFT_KIND, VDW_SWITCH_KIND;
};

#if ENSEMBLE == GEMC || ENSEMBLE == NPT

//Items that effect the system interactions and/or identity, e.g. Temp.
struct GEMCKind {
  uint kind;
  real pressure;
};

#endif


struct Step {
  ulong total, equil, adjustment, pressureCalcFreq;
  bool pressureCalc;
};

//Holds the percentage of each kind of move for this ensemble.
struct MovePercents {
  real displace, rotate, intraSwap, intraMemc, regrowth, crankShaft;
#ifdef VARIABLE_VOLUME
  real volume;
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  real transfer, memc;
#endif
};

struct ElectroStatic {
  bool readEwald;
  bool readElect;
  bool readCache;
  bool enable;
  bool ewald;
  bool cache;
  bool cutoffCoulombRead[BOX_TOTAL];
  real tolerance;
  real oneFourScale;
  real dielectric;
  real cutoffCoulomb[BOX_TOTAL];
  ElectroStatic(void)
  {
    std::fill_n(cutoffCoulombRead, BOX_TOTAL, false);
    std::fill_n(cutoffCoulomb, BOX_TOTAL, 0.0);
  }
};

struct Volume {
  bool hasVolume, cstArea, cstVolBox0;
  bool readCellBasis[BOX_TOTAL][3];
  XYZArray axis[BOX_TOTAL];
  Volume(void) : hasVolume(false), cstArea(false), cstVolBox0(false)
  {
    for(uint b = 0; b < BOX_TOTAL; ++b) {
      axis[b] = XYZArray(3);
      std::fill_n(readCellBasis[b], 3, false);
    }
  }

  bool ReadCellBasis() const
  {
    for(uint b = 0; b < BOX_TOTAL; ++b) {
      for(uint i = 0; i < 3; i++) {
        if(!readCellBasis[b][i])
          return false;
      }
    }
    return true;
  }
};

//If particle number varies (e.g. GCMC, GEMC) load in parameters for
//configurational bias
struct GrowNonbond {
  uint first, nth;
};

//If particle number varies (e.g. GCMC, GEMC) load in parameters for
//configurational bias
struct GrowBond {
  uint ang, dih;
};

struct CBMC {
  GrowNonbond nonbonded;
  GrowBond bonded;
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
  MEMCVal(void)
  {
    MEMC1 = MEMC2 = MEMC3 = false;
    readVol = readRatio = readSmallBB = false;
    readLargeBB = readSK = readLK = false;
  }
};

#if ENSEMBLE == GCMC
struct ChemicalPotential {
  bool isFugacity;
  std::map<std::string, real> cp;
};
#endif
struct SystemVals {
  ElectroStatic elect;
  Temperature T;
  FFValues ff;
  Exclude exclude;
  Step step;
  MovePercents moves;
  Volume volume; //May go unused
  CBMC cbmcTrials;
  MEMCVal memcVal, intraMemcVal;
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
  bool operator()(void)
  {
    return enable;
  }
};

struct UniqueStr { /* : ReadableBase*/
  std::string val;
};

struct HistFiles { /* : ReadableBase*/
  std::string histName, number, letter, sampleName;
  uint stepsPerHistSample;
};

//Files for output.
struct OutFiles {
  FileNames<BOX_TOTAL> pdb;
  FileName psf, seed;
  HistFiles hist;
};
struct Settings {
  EventSettings block, hist;
  UniqueStr uniqueStr;
};

//Enables for each variable that can be tracked
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
  SysState state, restart;
  Statistics statistics;
  EventSettings console, checkpoint;
};

}

class ConfigSetup
{
public:
  config_setup::Input in;
  config_setup::Output out;
  config_setup::SystemVals sys;
  ConfigSetup(void);
  void Init(const char *fileName);

private:
  void fillDefaults(void);
  bool checkBool(string str);
  bool CheckString(string str1, string str2);
  void verifyInputs(void);
  InputFileReader reader;

  //Names of config file.
  static const char defaultConfigFileName[]; // "in.dat"
  static const char configFileAlias[];       // "GO-MC Configuration File"
};

#endif
