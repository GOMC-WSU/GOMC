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
      ewaldDriven = statV.forcefield.ewald;
      originalWolfKind = statV.forcefield.GetWolfKind();
      originalCoulKind = statV.forcefield.GetCoulKind();
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            orignalWolfAlpha[b] = statV.forcefield.GetWolfAlpha(b);
            if(sysVals.wolfCal.wolfAlphaRangeRead[b]){
                  wolfAlphaStart[b] = sysVals.wolfCal.wolfAlphaStart[b];
                  wolfAlphaEnd[b] = sysVals.wolfCal.wolfAlphaEnd[b];
                  wolfAlphaDelta[b] = sysVals.wolfCal.wolfAlphaDelta[b];
                  alphaSize[b] = 0;
                  for (double a = wolfAlphaStart[b]; a <= wolfAlphaEnd[b]; a+=wolfAlphaDelta[b]){
                        alphaSize[b] += 1;
                  }
            }
      }
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            if(sysVals.wolfCal.wolfAlphaRangeRead[b]){
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              sumRelativeErrorVec[b][wolfKind][coulKind].resize(alphaSize[b]);
                              //relativeErrorVec[b][wolfKind][coulKind].resize(alphaSize[b]);
                              relativeError[b][wolfKind][coulKind] = new double[alphaSize[b]];
                        }
                  }
            }
      }
}

WolfCalibrationOutput::~WolfCalibrationOutput()
{
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        if (relativeError[b][wolfKind][coulKind] != NULL)
                              delete relativeError[b][wolfKind][coulKind];
                  }
            }
      }
}

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {

      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      if(enableOut) {
            WriteHeader();
            WriteGraceParFile();
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
                              /*
                              std::string firstRow = "";
                              firstRow += "Step#\t";
                              // We skip the reference r cut with reference alpha.
                              // r = 0, a = 0
                              // So there are no duplicate columns.
                              for (double a = wolfAlphaStart[b]; a <= wolfAlphaEnd[b]; a+=wolfAlphaDelta[b]){
                                    firstRow += std::to_string(a);
                                    firstRow += " ";
                              }
                              outF << firstRow;
                              outF << std::endl;
                              */
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
            outF.open(("WOLF_CALIBRATION_BOX_" + std::to_string(b) + ".par").c_str(), std::ofstream::out);
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

void WolfCalibrationOutput::DoOutput(const ulong step) {

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open(("WOLF_CALIBRATION_BOX_" + std::to_string(b) + ".dat").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  for (uint i = 0; i < alphaSize[b]; ++i) {
                        double a = wolfAlphaStart[b] + i*wolfAlphaDelta[b];
                        std::string firstRow = "";
                        firstRow += std::to_string(a) + "\t";
                        for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                              for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){ 
                                    double min_err = 100.00*((sumRelativeErrorVec[b][wolfKind][coulKind][i].mean()-ewaldAvg[b].mean())/ewaldAvg[b].mean());
                                    firstRow += std::to_string(min_err) + "\t";
                                    //firstRow += std::to_string(sumRelativeError[b][wolfKind][coulKind][i]/numSamples) + "\t";
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

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        outF.open((getFileName(b, wolfKind, coulKind, uniqueName)+".dat").c_str(), std::ofstream::app);            
                        if (outF.is_open()) {
                              std::string row = "";
                              row += GetString(step);
                              row += "\t";
                              for (uint i = 0; i < alphaSize[b]; ++i) {
                                    row += std::to_string(100.00*relativeError[b][wolfKind][coulKind][i]);
                                    row += "\t";
                              }
                              outF << row << std::endl;
                        } else {
                              std::cerr << "Unable to write to file \"" <<  name << "\" "
                                          << "(Wolf Calibration file)" << std::endl;
                        }
                        outF.close();
                  }
            }
      }
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            outF.open(("WOLF_CALIBRATION_BOX_" + std::to_string(b) + "_BEST_ALPHAS.csv").c_str(), std::ofstream::out);
            if (outF.is_open()) {
                  std::string firstRow = "";
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){ 
                              int min_i = 0;
                              double min_err = std::numeric_limits<double>::max();
                              for (uint i = 0; i < alphaSize[b]; ++i) {
                                    //printf("vec avg %f std %f z %f\n", sumRelativeErrorVec[b][wolfKind][coulKind][i].mean(), sumRelativeErrorVec[b][wolfKind][coulKind][i].sd(), std::abs(sumRelativeErrorVec[b][wolfKind][coulKind][i].mean()/sumRelativeErrorVec[b][wolfKind][coulKind][i].sd()));
                                    //double err = std::abs(sumRelativeErrorVec[b][wolfKind][coulKind][i].mean()/sumRelativeErrorVec[b][wolfKind][coulKind][i].sd());
                                    //double err = std::abs(sumRelativeErrorVec[b][wolfKind][coulKind][i].mean()-ewaldAvg[b].mean());
                                    double mean1 = ewaldAvg[b].mean();
                                    double sd1 = ewaldAvg[b].sd();
                                    double n = ewaldAvg[b].count();

                                    double mean2 = sumRelativeErrorVec[b][wolfKind][coulKind][i].mean();
                                    double sd2 = sumRelativeErrorVec[b][wolfKind][coulKind][i].sd();
                                    double m = sumRelativeErrorVec[b][wolfKind][coulKind][i].count();

                                    double t_test = (mean1 - mean2) / sqrt((sd1 * sd1) / n + (sd2 * sd2) / m);
                                    if (std::abs(t_test) < min_err){
                                          min_err = std::abs(t_test);
                                          min_i = i;
                                    }
                              }
                              double best_a = wolfAlphaStart[b] + min_i*wolfAlphaDelta[b];
                              std::string title = WOLF_KINDS[wolfKind] + " " + COUL_KINDS[coulKind];
                              firstRow += title + "\t" + std::to_string(best_a) + "\n";
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

}

/*
void WolfCalibrationOutput::DoOutput(const ulong step) {

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  statValRef.forcefield.SetWolfKind(wolfKind);
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        statValRef.forcefield.SetCoulKind(coulKind);
                        outF.open((getFileName(b, wolfKind, coulKind, uniqueName)+".dat").c_str(), std::ofstream::out);
                        if (outF.is_open()) {
                              for (uint i = 0; i < alphaSize[b]; ++i) {
                                    double a = wolfAlphaStart[b] + i*wolfAlphaDelta[b];
                                    std::string firstRow = "";
                                    firstRow += std::to_string(a) + "\t" + std::to_string(sumRelativeError[b][wolfKind][coulKind][i]/numSamples) + "\n";
                                    outF << firstRow;
                              }
                              outF << std::endl;
                        } else {
                              std::cerr << "Unable to write to file \"" <<  name << "\" "
                                          << "(Wolf Calibration file)" << std::endl;
                        }
                        outF.close();
                  }
            }
      }
}
*/
#include <float.h>

void WolfCalibrationOutput::Sample(const ulong step) {
      if (!((printOnFirstStep && step == startStep) ||
        (enableOut && ((step + 1) % stepsPerOut == 0) || forceOutput)))
            return;
      ++numSamples;


      SystemPotential ewaldRef;
      if (ewaldDriven) {
            ewaldRef = calcEn.SystemTotal();
            ewaldRef.Total();
            // Swap wolf and ewald
            std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);
      } else {
            //SystemPotential sanity = calcEn.SystemTotal();
            std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);
            #ifdef GOMC_CUDA
            statValRef.forcefield.particles->updateWolfEwald();
            #endif
            ewaldRef = calcEn.WolfCalSystemTotal();
            ewaldRef.Total();
            std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);

      }

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            //printf("EwAtStep %lu %.*e\n", step, Digs, ewaldRef.boxEnergy[b].totalElect);
            ewaldAvg[b].add_value(ewaldRef.boxEnergy[b].totalElect);
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  statValRef.forcefield.SetWolfKind(wolfKind);
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        statValRef.forcefield.SetCoulKind(coulKind);
                        for (uint i = 0; i < alphaSize[b]; ++i) {
                              double a = wolfAlphaStart[b] + i*wolfAlphaDelta[b];
                              // Wolf class has references to these forcefield values
                              statValRef.forcefield.SetWolfAlphaAndWolfFactors(a, b);
                              #ifdef GOMC_CUDA
                              statValRef.forcefield.particles->updateWolfEwald();
                              #endif
                              SystemPotential wolfTot;
                              if (ewaldDriven) {
                                    wolfTot = calcEn.WolfCalSystemTotal();
                              } else {
                                    wolfTot = calcEn.SystemTotal();
                              }
                              //printf("WoAtStep %lu %d %d %f %.*e\n", step, wolfKind, coulKind, a, Digs, wolfTot.boxEnergy[b].totalElect);
                              sumRelativeErrorVec[b][wolfKind][coulKind][i].add_value(wolfTot.boxEnergy[b].totalElect);
                              //relativeErrorVec[b][wolfKind][coulKind][i].push_back((wolfTot.boxEnergy[b].totalElect-ewaldRef.boxEnergy[b].totalElect)/ewaldRef.boxEnergy[b].totalElect);
                              relativeError[b][wolfKind][coulKind][i] = ((wolfTot.boxEnergy[b].totalElect-ewaldRef.boxEnergy[b].totalElect)/ewaldRef.boxEnergy[b].totalElect);
                        }
                  }
            }
      }
      if (ewaldDriven){
            std::swap(statValRef.forcefield.ewald, statValRef.forcefield.wolf);
            #ifdef GOMC_CUDA
            statValRef.forcefield.particles->updateWolfEwald();
            #endif
      } else {
            statValRef.forcefield.SetCoulKind(originalCoulKind);
            statValRef.forcefield.SetWolfKind(originalWolfKind);
            for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
                  statValRef.forcefield.SetWolfAlphaAndWolfFactors(orignalWolfAlpha[b], b);
            }
            #ifdef GOMC_CUDA
            statValRef.forcefield.particles->updateWolfEwald();
            #endif
      }
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