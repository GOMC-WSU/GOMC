/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "FreeEnergyOutput.h"

#include <iomanip>
#include <sstream>

#include "ConfigSetup.h"
#include "OutConst.h"
#include "PDBConst.h"

FreeEnergyOutput::FreeEnergyOutput(OutputVars &v, System &sys)
    : calcEn(sys.calcEnergy), freeEnVal(sys.statV.freeEnVal),
      lambdaRef(sys.lambdaRef) {
  this->var = &v;
  for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
    energyDiff[b] = NULL;
  }
#if ENSEMBLE == NPT
  // unit is K * molecule / A3
  imposedP = sys.statV.pressure;
#endif
}

void FreeEnergyOutput::Init(pdb_setup::Atoms const &atoms,
                            config_setup::Output const &output) {
  stepsPerSample = freeEnVal.frequency;
  stepsPerOut = freeEnVal.frequency;
  enableOut = freeEnVal.enable;
  lambdaSize = freeEnVal.lambdaVDW.size();
  iState = freeEnVal.iState;

  if (enableOut) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      std::stringstream sstrm;
      std::string strKind, fileName;
      sstrm << (b);
      sstrm >> strKind;
      fileName = "Free_Energy_BOX_";
      fileName += strKind;
      fileName += "_";
      fileName += uniqueName;
      fileName += ".dat";
#if GOMC_LIB_MPI
      name[b] = pathToReplicaOutputDirectory + fileName;
#else
      name[b] = fileName;
#endif
      outF[b].open(name[b].c_str(), std::ofstream::out);
      energyDiff[b] = new Energy[lambdaSize];
    }
    WriteHeader();
  }
}

std::string FreeEnergyOutput::GetString(double a, uint p) {
  std::stringstream sstrm;
  std::string tempStr;
  sstrm << std::fixed << std::setprecision(p) << (a);
  // sstrm.precision(p);
  // sstrm >> tempStr;
  tempStr = sstrm.str();
  return tempStr;
}

FreeEnergyOutput::~FreeEnergyOutput() {
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) {
      outF[b].close();
    }

    if (energyDiff[b] != NULL) {
      delete[] energyDiff[b];
    }
  }
}

void FreeEnergyOutput::Sample(const ulong step) {
  if ((step + 1) % stepsPerSample == 0) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      CalculateFreeEnergy(b);
    }
  }
}

void FreeEnergyOutput::DoOutput(const ulong step) {
  // Write to histogram file, We don't check the equilibrium.
  GOMC_EVENT_START(1, GomcProfileEvent::FREE_ENERGY_OUTPUT);
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) {
      PrintData(b, step + 1);
    } else {
      std::cerr << "Unable to write to file \"" << name[b] << "\" "
                << "(Free Energy file)" << std::endl;
      outF[b].close();
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::FREE_ENERGY_OUTPUT);
}

void FreeEnergyOutput::DoOutputRestart(const ulong step) {}

void FreeEnergyOutput::PrintData(const uint b, const ulong step) {
  outF[b] << std::setw(11) << std::left << step << " ";
  outF[b] << std::setw(25) << std::right << std::fixed << Etotal << " ";
  outF[b] << std::setw(25) << dUdL_Coulomb[b].Total() << " ";
  outF[b] << std::setw(25) << dUdL_VDW[b].Total() << " ";

  for (uint i = 0; i < lambdaSize; i++) {
    outF[b] << std::setw(25) << energyDiff[b][i].Total() << " ";
  }
#if ENSEMBLE == NVT
  if (var->pressureCalc) {
    outF[b] << std::setw(25) << PV;
  }
#elif ENSEMBLE == NPT
  outF[b] << std::setw(25) << PV;
#endif
  outF[b] << std::endl;
}

void FreeEnergyOutput::WriteHeader(void) {
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) {
      std::string toPrint = "";
      toPrint += "#T = ";
      toPrint += GetString(var->T_in_K, 4);
      toPrint += "(K), Lambda State ";
      toPrint += GetString(iState, 0);
      toPrint += ": (lambda Coulomb, lambda VDW) = (";
      toPrint += GetString(freeEnVal.lambdaCoulomb[iState], 4);
      toPrint += ",";
      toPrint += GetString(freeEnVal.lambdaVDW[iState], 4);
      toPrint += ")\n";
      outF[b] << toPrint;

      // We care about format
      outF[b] << std::setw(11) << std::left << "#Steps"
              << " ";
      outF[b] << std::setw(25) << std::right << "Total_En(kJ/mol)"
              << " ";
      toPrint = "dU/dL(Coulomb=";
      toPrint += GetString(freeEnVal.lambdaCoulomb[iState], 4);
      toPrint += ")";
      outF[b] << std::setw(25) << std::right << toPrint << " ";
      toPrint = "dU/dL(VDW=";
      toPrint += GetString(freeEnVal.lambdaVDW[iState], 4);
      toPrint += ")";
      outF[b] << std::setw(25) << std::right << toPrint << " ";

      std::string fixStr = "DelE(L->(";
      for (uint i = 0; i < lambdaSize; i++) {
        toPrint = fixStr;
        toPrint += GetString(freeEnVal.lambdaCoulomb[i], 4);
        toPrint += ",";
        toPrint += GetString(freeEnVal.lambdaVDW[i], 4);
        toPrint += "))";
        outF[b] << std::setw(25) << std::right << toPrint << " ";
      }
#if ENSEMBLE == NVT
      if (var->pressureCalc) {
        outF[b] << std::setw(25) << std::right << "PV(kJ/mol)";
      }
#elif ENSEMBLE == NPT
      outF[b] << std::setw(25) << std::right << "PV(kJ/mol)";
#endif
      outF[b] << std::endl;
      outF[b] << std::setprecision(10);
      outF[b].setf(std::ios_base::right, std::ios_base::adjustfield);
    } else
      std::cerr << "Unable to write to file \"" << name[b] << "\" "
                << "(Free Energy file)" << std::endl;
  }
}

void FreeEnergyOutput::CalculateFreeEnergy(const uint b) {
#if ENSEMBLE == NVT
  if (var->pressureCalc) {
    PV = var->pressure[b] * var->volumeRef[b] * unit::BAR_TO_K_MOLECULE_PER_A3;
    PV *= unit::K_TO_KJ_PER_MOL;
  }
#elif ENSEMBLE == NPT
  // no need to convert pressure (bar) to K * molecule /A3
  PV = imposedP * var->volumeRef[b] * unit::K_TO_KJ_PER_MOL;
#endif
  Etotal = var->energyRef[b].Total() * unit::K_TO_KJ_PER_MOL;
  uint molIndex = lambdaRef.GetMolIndex(b);
  // Reset the energy value
  dUdL_VDW[b].Zero();
  dUdL_Coulomb[b].Zero();
  for (uint i = 0; i < lambdaSize; i++) {
    energyDiff[b][i].Zero();
  }
  // Calculate delta E and dE/dLambda for all lambda states
  calcEn.EnergyChange(energyDiff[b], dUdL_VDW[b], dUdL_Coulomb[b],
                      freeEnVal.lambdaVDW, freeEnVal.lambdaCoulomb, iState,
                      molIndex, b);
  // Convert to kJ/mol
  dUdL_VDW[b] *= unit::K_TO_KJ_PER_MOL;
  dUdL_Coulomb[b] *= unit::K_TO_KJ_PER_MOL;
  for (uint i = 0; i < lambdaSize; i++) {
    energyDiff[b][i] *= unit::K_TO_KJ_PER_MOL;
  }
}
