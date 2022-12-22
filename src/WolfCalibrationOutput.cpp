/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "WolfCalibrationOutput.h"
#include "GOMC_Config.h"


WolfCalibrationOutput::WolfCalibrationOutput(System & sys, StaticVals & statV, config_setup::SystemVals const &sysVals):
sysRef(sys), calcEn(sys.calcEnergy), statValRef(statV)
{
      // This is neccessary to check for correctness of single point energy calculations.
      printOnFirstStep = true;
      numSamples = 0;
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            wolfAlphaStart[b] = sysVals.wolfCal.wolfAlphaStart[b];
            wolfAlphaEnd[b] = sysVals.wolfCal.wolfAlphaEnd[b];
            wolfAlphaDelta[b] = sysVals.wolfCal.wolfAlphaDelta[b];
            alphaSize[b] = 0;
            for (double a = wolfAlphaStart[b]; a <= wolfAlphaEnd[b]; a+=wolfAlphaDelta[b]){
                  alphaSize[b] += 1;
            }
      }
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        sumRelativeError[b][wolfKind][coulKind]  = new double[alphaSize[b]];
                        for (uint i = 0; i < alphaSize[b]; ++i) {
                              sumRelativeError[b][wolfKind][coulKind][i] = 0.0;
                        }
                  }
            }
      }
}

WolfCalibrationOutput::~WolfCalibrationOutput()
{
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) 
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind)
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind)
                        delete sumRelativeError[b][wolfKind][coulKind];
      
}

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {

      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      if(enableOut) {
            //WriteHeader();
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
                        } else {
                              std::cerr << "Unable to write to file \"" <<  name << "\" "
                                          << "(Wolf Calibration file)" << std::endl;
                        }
                        outF.close();
                  }
            }
      }
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
                                    firstRow += std::to_string(sumRelativeError[b][wolfKind][coulKind][i]/numSamples) + "\t";
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
 
void WolfCalibrationOutput::Sample(const ulong step) {
      if (!((printOnFirstStep && step == startStep) ||
        (enableOut && ((step + 1) % stepsPerOut == 0) || forceOutput)))
            return;
      ++numSamples;
      SystemPotential ewaldRef = calcEn.SystemTotal();
      ewaldRef.Total();

      // Swap wolf and ewald
      sysRef.SwapWolfAndEwaldPointers();
      bool tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;

      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
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
                              SystemPotential wolfTot = calcEn.SystemTotal();
                              sumRelativeError[b][wolfKind][coulKind][i] += (wolfTot.boxEnergy[b].totalElect-ewaldRef.boxEnergy[b].totalElect)/ewaldRef.boxEnergy[b].totalElect;
                        }
                  }
            }
      }
      sysRef.SwapWolfAndEwaldPointers();
      tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;
      #ifdef GOMC_CUDA
      statValRef.forcefield.particles->updateWolfEwald();
      #endif
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