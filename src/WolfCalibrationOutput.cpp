/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "WolfCalibrationOutput.h"
#include "GOMC_Config.h"
#include <limits>

WolfCalibrationOutput::WolfCalibrationOutput(System & sys, StaticVals & statV, config_setup::SystemVals const &sysVals):
sysRef(sys), calcEn(sys.calcEnergy), statValRef(statV)
{
      // This is neccessary to check for correctness of single point energy calculations.
      printOnFirstStep = true;
      //printOnFirstStep = false;
      numSamples = 0;
      for (int i = 0; i < BOXES_WITH_U_NB; ++i){
            wolfAlphaRangeRead[i]=sysVals.wolfCal.wolfAlphaRangeRead[i];
            wolfCutoffCoulombRangeRead[i]=sysVals.wolfCal.wolfCutoffCoulombRangeRead[i];
      }
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            if(wolfAlphaRangeRead[b]){
                  wolfAlphaStart[b] = sysVals.wolfCal.wolfAlphaStart[b];
                  wolfAlphaEnd[b] = sysVals.wolfCal.wolfAlphaEnd[b];
                  wolfAlphaDelta[b] = sysVals.wolfCal.wolfAlphaDelta[b];
                  alphaSize[b] = 0;
                  for (double a = wolfAlphaStart[b]; a <= wolfAlphaEnd[b]; a+=wolfAlphaDelta[b]){
                        alphaSize[b] += 1;
                  }
            }
            originalCutoffCoulomb[b] = statV.forcefield.rCutCoulomb[b];
            originalCutoffCoulombSq[b] = statV.forcefield.rCutCoulombSq[b];
            if(wolfCutoffCoulombRangeRead[b]){
                  wolfCutoffCoulombStart[b] = sysVals.wolfCal.wolfCutoffCoulombStart[b];
                  wolfCutoffCoulombEnd[b] = sysVals.wolfCal.wolfCutoffCoulombEnd[b];
                  wolfCutoffCoulombDelta[b] = sysVals.wolfCal.wolfCutoffCoulombDelta[b];
                  cutoffCoulombSize[b] = 0;
                  for (double a = wolfCutoffCoulombStart[b]; a <= wolfCutoffCoulombEnd[b]; a+=wolfCutoffCoulombDelta[b]){
                        cutoffCoulombSize[b] += 1;
                  }
            }
      }
}

WolfCalibrationOutput::~WolfCalibrationOutput()
{

}

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {

      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      if(enableOut) {
            //WriteHeader();
            //WriteGraceParFile();
            for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
                  if(wolfAlphaRangeRead[b] && wolfCutoffCoulombRangeRead[b]){
                        for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                              for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                                    sumRelativeErrorVec[b][wolfKind][coulKind].resize(alphaSize[b]*cutoffCoulombSize[b]);
                              }
                        }
                  }
            }
      }
}

std::string WolfCalibrationOutput::getFileName(int b, int wolfKind, int coulKind, std::string uniqueName){
      std::stringstream sstrm;
      std::string strKind, fileName, name;
      sstrm << (b);
      sstrm >> strKind;
      fileName = "Wolf_Calibration_";
      fileName += WOLF_KINDS[wolfKind];
      fileName += "_";
      fileName += COUL_KINDS[coulKind];
      fileName += "_BOX_";
      fileName += strKind;
      fileName += "_";
      fileName += uniqueName;
      #if GOMC_LIB_MPI
            name = pathToReplicaOutputDirectory + fileName;
      #else
            name = fileName;
      #endif
      return name;
}

void WolfCalibrationOutput::WriteHeader()
{
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        outF.open((getFileName(b, wolfKind, coulKind, uniqueName)+".dat").c_str(), std::ofstream::out);
                        if (outF.is_open()) {

                        } else {
                              std::cerr << "Unable to write to file \"" <<  name << "\" "
                                          << "(Wolf Calibration file)" << std::endl;
                        }
                        outF.close();
                  }
            }
      }
}

template<typename T>
double getAverage(std::vector<T> const& v) {
      if (v.empty()) {
            return 0;
      }
      return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<typename T>
double getSTD(std::vector<T> const& v) {
      if (v.empty()) {
            return 0;
      }
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      double mean = sum / v.size();

      std::vector<double> diff(v.size());
      std::transform(v.begin(), v.end(), diff.begin(),
                  std::bind2nd(std::minus<double>(), mean));
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      double stdev = std::sqrt(sq_sum / v.size());
      return stdev;
}

void WolfCalibrationOutput::WriteGraceParFile()
{
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open((uniqueName + "_WOLF_CALIBRATION_BOX_" + std::to_string(b) + ".par").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  int counter = 0;
                  std::string firstRow = "";
                  firstRow += "title \"Comparing Wolf Models\"\n";
                  firstRow += "xaxis label \"Alpha\"\n";
                  firstRow += "yaxis label \"Relative Error\"\n";
                  firstRow += "TITLE SIZE 2 \n";
                  firstRow += "LEGEND .8,.45\n";
                  firstRow += "with g0\n";
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              std::string title = WOLF_KINDS[wolfKind] + " " + COUL_KINDS[coulKind];
                              firstRow += "\ts";
                              firstRow += GetString(counter);
                              firstRow += " legend \"" + title + "\"\n";
                              counter += 1;
                        }
                  }
                  outF << firstRow;
                  outF << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
            outF.close();
      }
}


void WolfCalibrationOutput::WriteGraceParFileWRcut()
{
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open((uniqueName + "_WOLF_CALIBRATION_BOX_" + std::to_string(b) + ".par").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  int counter = 0;
                  std::string firstRow = "";
                  firstRow += "title \"Comparing Wolf Models\"\n";
                  firstRow += "xaxis label \"Alpha\"\n";
                  firstRow += "yaxis label \"Relative Error\"\n";
                  firstRow += "TITLE SIZE 2 \n";
                  firstRow += "LEGEND .8,.45\n";
                  firstRow += "with g0\n";
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              int min_rcut_index=mapWK_CK_BOX_to_bestRCutIndex[std::make_tuple(wolfKind, coulKind, b)];
                              double best_rcut = wolfCutoffCoulombStart[b] + min_rcut_index*wolfCutoffCoulombDelta[b];
                              std::ostringstream oss;
                              oss << std::fixed << std::setprecision(3) << best_rcut;
                              std::string result = oss.str();
                              std::string title = WOLF_KINDS[wolfKind] + ", " + COUL_KINDS[coulKind] + ", RCut " + result;
                              // Insert values into the map
                              firstRow += "\ts";
                              firstRow += GetString(counter);
                              firstRow += " legend \"" + title + "\"\n";
                              counter += 1;
                        }
                  }
                  outF << firstRow;
                  outF << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
            outF.close();
      }
}

void WolfCalibrationOutput::DoOutput(const ulong step) {
      //if (((step+1) < stepsTillEquil) || !(enableOut && ((step + 1) % stepsPerOut == 0)))
      //      return;
      CalculateGrid();
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        outF.open((getFileName(b, wolfKind, coulKind, uniqueName)+".dat").c_str(), std::ofstream::out);            
                        if (outF.is_open()) {
                              std::string firstRow = "";
                              for (uint RCutIndex = 0; RCutIndex < cutoffCoulombSize[b]; ++RCutIndex) {
                                    firstRow += ",";
                                    double rCutCoulomb = wolfCutoffCoulombStart[b] + RCutIndex*wolfCutoffCoulombDelta[b];
                                    firstRow += std::to_string(rCutCoulomb);
                              }
                              outF << firstRow << std::endl;
                              for (uint alphaIndex = 0; alphaIndex < alphaSize[b]; ++alphaIndex) {
                                    std::string dataRow = "";
                                    double alpha = wolfAlphaStart[b] + alphaIndex*wolfAlphaDelta[b];
                                    dataRow += std::to_string(alpha);
                                    for (uint RCutIndex = 0; RCutIndex < cutoffCoulombSize[b]; ++RCutIndex) {
                                          dataRow += ",";
                                          double relativeError = 100.00*((sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].mean()-ewaldAvg[b].mean())/ewaldAvg[b].mean());
                                          dataRow += std::to_string(relativeError);
                                    }
                                    outF << dataRow << std::endl;
                              }
                              outF.close();
                        }
                  }
            }
      }

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open((uniqueName + "_WOLF_CALIBRATION_BOX_" + std::to_string(b) + "_BEST_ALPHAS.csv").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  std::string firstRow = "WolfKind,CoulKind,Rcut,Alpha,Error\n";
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){ 
                              int min_i = 0;
                              int min_a_index = 0;
                              int min_rcut_index = 0;
                              double min_err = 0.0;
                              double min_t_stat = std::numeric_limits<double>::max();
                              for (uint alphaIndex = 0; alphaIndex < alphaSize[b]; ++alphaIndex) {
                                    double alpha = wolfAlphaStart[b] + alphaIndex*wolfAlphaDelta[b];
                                    for (uint RCutIndex = 0; RCutIndex < cutoffCoulombSize[b]; ++RCutIndex) {
                                          double rCutCoulomb = wolfCutoffCoulombStart[b] + RCutIndex*wolfCutoffCoulombDelta[b];
                                          double err = std::abs(sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].mean()-ewaldAvg[b].mean());
                                          double mean1 = ewaldAvg[b].mean();
                                          double sd1 = ewaldAvg[b].sd();
                                          double n = ewaldAvg[b].count();

                                          double mean2 = sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].mean();
                                          double sd2 = sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].sd();
                                          double m = sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].count();

                                          double t_test;
                                          if (n == 1)
                                                t_test = (mean1 - mean2);
                                          else
                                                t_test = (mean1 - mean2) / sqrt((sd1 * sd1) / n + (sd2 * sd2) / m);

                                          if (std::abs(t_test) < min_t_stat){
                                                min_t_stat = std::abs(t_test);
                                                min_err = err;
                                                min_a_index = alphaIndex;
                                                min_rcut_index = RCutIndex;
                                          }
                                    }
                              }
                              double best_a = wolfAlphaStart[b] + min_a_index*wolfAlphaDelta[b];
                              double best_rcut = wolfCutoffCoulombStart[b] + min_rcut_index*wolfCutoffCoulombDelta[b];
                              std::string title = WOLF_KINDS[wolfKind] + "," + COUL_KINDS[coulKind];
                              firstRow += title + "," + std::to_string(best_rcut) + "," + std::to_string(best_a) + "," + std::to_string(min_err) + "\n";
                              // Insert values into the map
                              mapWK_CK_BOX_to_bestRCutIndex[std::make_tuple(wolfKind, coulKind, b)] = min_rcut_index;
                        }
                  }
                  firstRow += "\n";
                  outF << firstRow;
                  outF << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
            outF.close();
      }

      
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open((uniqueName + "_WOLF_CALIBRATION_BOX_" + std::to_string(b) + ".dat").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  for (uint alphaIndex = 0; alphaIndex < alphaSize[b]; ++alphaIndex) {
                        double a = wolfAlphaStart[b] + alphaIndex*wolfAlphaDelta[b];
                        std::string firstRow = "";
                        firstRow += std::to_string(a) + "\t";
                        for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                              for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                                    int min_rcut_index = mapWK_CK_BOX_to_bestRCutIndex[std::make_tuple(wolfKind, coulKind, b)];
                                    double min_err = 100.00*((sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(min_rcut_index, alphaIndex, b)].mean()-ewaldAvg[b].mean())/ewaldAvg[b].mean());
                                    firstRow += std::to_string(min_err) + "\t";
                                    //firstRow += std::to_string(sumRelativeError[b][i]/numSamples) + "\t";
                              }
                        }
                        firstRow += "\n";
                        outF << firstRow;
                  }
                  outF << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
            outF.close();
      }     
      WriteGraceParFileWRcut();
}

#include <float.h>

void WolfCalibrationOutput::CalculateGrid() {
      ++numSamples;
      SystemPotential ewaldRef;
      ewaldRef = calcEn.SystemTotal();
      ewaldRef.Total();
      // Swap wolf and ewald
      std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);
      std::swap(sysRef.calcEwald, sysRef.calcWolf);
      calcEn.UpdateEwald();

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            //printf("EwAtStep %lu %.*e\n", step, Digs, ewaldRef.boxEnergy[b].totalElect);
            ewaldAvg[b].add_value(ewaldRef.boxEnergy[b].totalElect);
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        //statValRef.forcefield.dsf = coulKind;
                        statValRef.forcefield.SetWolfMethod(wolfKind,coulKind);
                        for (uint alphaIndex = 0; alphaIndex < alphaSize[b]; ++alphaIndex) {
                              double alpha = wolfAlphaStart[b] + alphaIndex*wolfAlphaDelta[b];
                              for (uint RCutIndex = 0; RCutIndex < cutoffCoulombSize[b]; ++RCutIndex) {
                                    double rCutCoulomb = wolfCutoffCoulombStart[b] + RCutIndex*wolfCutoffCoulombDelta[b];
                                    // Wolf class has references to these forcefield values
                                    statValRef.forcefield.SetWolfAlphaAndWolfFactors(rCutCoulomb, alpha, b);
                                    //#ifdef GOMC_CUDA
                                    //statValRef.forcefield.particles->updateWolfEwald();
                                    //#endif
                                    SystemPotential wolfTot = calcEn.SystemTotal();
                                    wolfTot.Total();
                                    sumRelativeErrorVec[b][wolfKind][coulKind][GetIndex(RCutIndex, alphaIndex, b)].add_value(wolfTot.boxEnergy[b].totalElect);
                              }
                        }
                  }
            }
      }
      // Swap wolf and ewald
      std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);
      std::swap(sysRef.calcEwald, sysRef.calcWolf);
      calcEn.UpdateEwald();
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            statValRef.forcefield.SetRCutCoulomb(originalCutoffCoulomb[b], b);
      }
      //#ifdef GOMC_CUDA
      //statValRef.forcefield.particles->updateWolfEwald();
      //#endif

}


std::string WolfCalibrationOutput::GetString(double a, uint p)
{
      std::stringstream sstrm;
      std::string tempStr;
      sstrm << std::fixed << std::setprecision(p) << (a);
      //sstrm.precision(p);
      //sstrm >> tempStr;
      tempStr = sstrm.str();
      return tempStr;
}

std::string WolfCalibrationOutput::GetString(ulong step)
{
      std::stringstream ss;
      ss << step;
      return ss.str();
}

// One row has a constant alpha and a varying Rcut
int WolfCalibrationOutput::GetIndex(int RCutIndex, int alphaIndex, int b){
      return alphaIndex*cutoffCoulombSize[b] + RCutIndex;
}
