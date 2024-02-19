/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef WOLF_CALIBRATION_OUTPUT_H
#define WOLF_CALIBRATION_OUTPUT_H

#include "OutputAbstracts.h"
#include <iostream>
#include "GOMC_Config.h"
#include "System.h"
#include "CalculateEnergy.h"
#include "Welford.h"
#include <map>
#include <tuple>


#include <string.h>
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
const std::string COUL_DSP_STRING = "DSP";
const std::string COUL_DSF_STRING = "DSF";
const std::string WOLF_RAHBARI_STRING = "RAHBARI";
const std::string WOLF_WAIBEL2018_STRING = "WAIBEL2018";
const std::string WOLF_WAIBEL2019_STRING = "WAIBEL2019";

enum COUL_KINDS_ENUMS { COUL_DSP, COUL_DSF };
enum WOLF_KINDS_ENUMS { WOLF_RAHBARI, WOLF_WAIBEL2018, WOLF_WAIBEL2019 };

const std::string COUL_KINDS[2] = {COUL_DSP_STRING, COUL_DSF_STRING};
const std::string WOLF_KINDS[3] = {WOLF_RAHBARI_STRING, WOLF_WAIBEL2018_STRING, WOLF_WAIBEL2019_STRING};
const int COUL_TOTAL_KINDS = 2;
const int WOLF_TOTAL_KINDS = 3;


class WolfCalibrationOutput : public OutputableBase
{
public:
  WolfCalibrationOutput(System & sys, StaticVals & statV, config_setup::SystemVals const &sysVals);

  ~WolfCalibrationOutput();

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step) {}  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step);

private:
  std::string getFileName(int b, int wolfKind, int coulKind, std::string uniqueName);

  void WriteHeader();
  void WriteGraceParFile();
  void WriteGraceParFileWRcut();

  std::string GetString(double a, uint p);
  std::string GetString(ulong step);
  int GetIndex(int RCutIndex, int alphaIndex, int b);


  System & sysRef;
  StaticVals & statValRef;
  CalculateEnergy & calcEn;
  uint stepsPerSample;


  std::ofstream outF;
  std::string name;
  
  double wolfAlphaStart[BOX_TOTAL];
  double wolfAlphaEnd[BOX_TOTAL];
  double wolfAlphaDelta[BOX_TOTAL];
  double wolfCutoffCoulombStart[BOX_TOTAL];
  double wolfCutoffCoulombEnd[BOX_TOTAL];
  double wolfCutoffCoulombDelta[BOX_TOTAL];
  int alphaSize[BOX_TOTAL];
  int cutoffCoulombSize[BOX_TOTAL];
  bool wolfAlphaRangeRead[2];
  bool wolfCutoffCoulombRangeRead[2];

  std::vector<Welford<double>> sumRelativeErrorVec[BOX_TOTAL];
  //std::vector<std::vector<double>> relativeErrorVec[BOX_TOTAL][WOLF_TOTAL_KINDS][COUL_TOTAL_KINDS];
  double *relativeError[BOX_TOTAL];
  std::map<std::tuple<int, int, int>, int> mapWK_CK_BOX_to_bestRCutIndex;

  Welford<double> ewaldAvg[BOX_TOTAL];
  int numSamples;
  double originalCutoffCoulomb[BOX_TOTAL];
  double originalCutoffCoulombSq[BOX_TOTAL];

};

#endif

