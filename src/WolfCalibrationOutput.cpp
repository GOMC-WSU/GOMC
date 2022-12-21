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
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
            wolfAlphaStart[b] = sysVals.wolfCal.wolfAlphaStart[b];
            wolfAlphaEnd[b] = sysVals.wolfCal.wolfAlphaEnd[b];
            wolfAlphaDelta[b] = sysVals.wolfCal.wolfAlphaDelta[b];
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
            for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              WriteHeader(b, wolfKind, coulKind);
                              WriteGraceParFile(b, wolfKind, coulKind);
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

void WolfCalibrationOutput::WriteHeader(uint b, uint wolfKind, uint coulKind)
{
      outF.open(getFileName(b, wolfKind, coulKind, uniqueName).c_str(), std::ofstream::out);
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


void WolfCalibrationOutput::WriteGraceParFile(uint b, uint wolfKind, uint coulKind)
{
      /*
      if (outFPar.is_open()) {
            int counter = 0;
            std::string firstRow = "";
            firstRow += "with g0\n";
            // We skip the reference r cut with reference alpha.
            // r = 0, a = 0
            // So there are no duplicate columns.
            for (int r = 0; r < wolfCalRef.numberOfRCuts[b]; ++r){
                  for (int a = 0; a < wolfCalRef.numberOfAlphas[b]; ++a){
                        firstRow += "\ts";
                        firstRow += GetString(counter);
                        firstRow += " legend \"(";
                        firstRow += GetString(wolfCalRef.rCutCoulomb[wolfCalRef.startOfNumRCuts[b]+r], 4);
                        firstRow += ", ";
                        firstRow += GetString(wolfCalRef.wolfAlpha[wolfCalRef.startOfNumAlphas[b]+a], 4);
                        firstRow += ")\"\n";
                        ++counter;
                  }
            }
            outFPar << firstRow;
            outFPar << std::endl;
      } else {
            std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                        << "(Wolf Calibration file)" << std::endl;
      }
      */
}

void WolfCalibrationOutput::DoOutput(const ulong step) {
      // Wolf driven GOMC cal
      // Buggy for NPT
      /*
      uint wolfKindOrig = sysRef.calcEwald->GetWolfKind();
      uint coulKindOrig = sysRef.calcEwald->GetCoulKind();
      calcEn.WolfCalibrationEnergy(&electrostaticEnergies[0]);

      // Calc reference epot
      sysRef.SwapWolfAndEwaldPointers();
      
      bool tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;
      //statValRef.forcefield.ewald != statValRef.forcefield.ewald;
      //statValRef.forcefield.wolf != statValRef.forcefield.wolf;
      
      // Set isVlugtWolf == false
      // Since ewald and wolf share a ff object
      statValRef.forcefield.SetWolfKind(0);
      // Ewald K Vectors may be not up to date with box dims
      // Since Wolf drives the simulation.
      //calcEn.ReInitKVectors();
      SystemPotential ewaldRef = calcEn.SystemTotal();
      ewaldRef.Total();
      sysRef.SwapWolfAndEwaldPointers();

      // Restore original wolf settings
      tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;

      // Restore inter wolf settings
      statValRef.forcefield.SetWolfKind(wolfKindOrig);
      statValRef.forcefield.SetCoulKind(coulKindOrig);

      // Restore intra wolf settings
      sysRef.calcEwald->SetWolfKind(wolfKindOrig);
      sysRef.calcEwald->SetCoulKind(coulKindOrig);

      */

      // Ewald driven GOMC cal
      /*
      SystemPotential ewaldRef = calcEn.SystemTotal();
      ewaldRef.Total();

      // Calc wolfcal epots
      sysRef.SwapWolfAndEwaldPointers();
      bool tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;

      calcEn.WolfCalibrationEnergy(&electrostaticEnergies[0]);

      sysRef.SwapWolfAndEwaldPointers();
      tmp = statValRef.forcefield.ewald;
      statValRef.forcefield.ewald = statValRef.forcefield.wolf;
      statValRef.forcefield.wolf = tmp;
    
      std::string row = "";
      row += GetString(step);
      row += "\t";
      
      for (uint box = 0; box < BOXES_WITH_U_NB; ++box) {       
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){  
                        std::string row = "";
                        row += GetString(step);
                        row += "\t";
                        // We skip the reference r cut with reference alpha.
                        // r = 0, a = 0
                        // So there are no duplicate columns.
                        for (int r = 0; r < wolfCalRef.numberOfRCuts[box]; ++r){
                              for (int a = 0; a < wolfCalRef.numberOfAlphas[box]; ++a){
                                    // If you dont use std::abs, double is converted to int 
                                    row += GetString((std::abs(ewaldRef.boxEnergy[box].TotalElect()) -  std::abs(electrostaticEnergies[wolfCalRef.GetIndex(box, wolfKind, coulKind, r, a)]))/ std::abs(ewaldRef.boxEnergy[box].TotalElect()), 16);
                                    row += "\t";
                              }
                        }
                        outF[box][wolfKind][coulKind] << row;
                        outF[box][wolfKind][coulKind] << std::endl;
                  }
            }
      }
      */
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