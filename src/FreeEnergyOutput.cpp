/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FreeEnergyOutput.h"
#include "PDBConst.h"
#include "OutConst.h"
#include "ConfigSetup.h"

#include <sstream>

FreeEnergyOutput::FreeEnergyOutput(OutputVars & v, System & sys) : calcEn(sys.calcEnergy), freeEnVal(sys.statV.freeEnVal), lambdaRef(sys.lambdaRef)
{
  this->var = &v;
  for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
    energyDiff[b] = NULL;
  }
}

void FreeEnergyOutput::Init(pdb_setup::Atoms const& atoms,
                     config_setup::Output const& output)
{
  stepsPerSample = freeEnVal.frequency;
  stepsPerOut = freeEnVal.frequency;
  enableOut = freeEnVal.enable;
  lambdaSize = freeEnVal.lambdaVDW.size();
  iState = freeEnVal.iState;

  if(enableOut) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      std::stringstream sstrm;  
      std::string strKind, fileName;
      sstrm << (b);
      sstrm >> strKind;
      fileName = "Free_Energy_";
      fileName += uniqueName;
      fileName += "_";
      fileName += strKind;
      fileName += ".dat";
      name[b] = fileName;
      outF[b].open(name[b].c_str(), std::ofstream::out);
      energyDiff[b] = new Energy[lambdaSize];
    }
    WriteHeader(); 
  }
}

FreeEnergyOutput::~FreeEnergyOutput()
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if(outF[b].is_open()) {
      outF[b].close();
    }

    if (energyDiff[b] != NULL) {
      delete[] energyDiff[b];
    }
  }
}

void FreeEnergyOutput::Sample(const ulong step)
{
  if ((step + 1) % stepsPerSample == 0) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      CalculateFreeEnergy(b);
    }
  }
}

void FreeEnergyOutput::DoOutput(const ulong step)
{
  //Write to histogram file, We dont check the equilibrium.
  if ((step + 1) % stepsPerOut == 0) {
    for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
      if (outF[b].is_open()) {
          PrintData(b, step + 1);
      } else {
        std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                  << "(Free Energy file)" << std::endl;
        outF[b].close();
      }
    }
  }
}

void FreeEnergyOutput::PrintData(const uint b, const uint step)
{
  outF[b] << std::setw(11) << std::left << step << " ";
  outF[b] << std::setw(25) << std::right << std::fixed << Etotal << " ";
  outF[b] << std::setw(25) << dUdL_Coulomb[b].Total();
  outF[b] << std::setw(25) << dUdL_VDW[b].Total() ;

  for(uint i = 0; i < lambdaSize; i++) {
    outF[b] << std::setw(25) << energyDiff[b][i].Total();
  }
  outF[b] << std::setw(25) << PV << std::endl;
}

void FreeEnergyOutput::WriteHeader(void)
{
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if (outF[b].is_open()) {
      outF[b] << "#T = " << std::setprecision(4) <<  std::fixed << 
                var->T_in_K << "(K), " << "Lambda State " <<
                iState << ": (lambda Coulomb, lambda VDW) = (" << 
                freeEnVal.lambdaCoulomb[iState] << "," << 
                freeEnVal.lambdaVDW[iState] << ")\n";

      //We care about format
      outF[b] << std::setw(11) << std::left << "#Steps" << " ";
      outF[b] << std::setw(25) << std::right << "Total_En(kJ/mol)" << " ";
      outF[b] << std::setw(25) << std::right << 
                "dU/dL(Coulomb = " << std::setprecision(4) <<
                freeEnVal.lambdaCoulomb[iState] << ") " ;
      outF[b] << std::setw(25) << "dU/dL(VDW=" << std::setprecision(4) <<
                freeEnVal.lambdaVDW[iState] << ") ";
      std::string fixStr = "DeltaE (L to (";
      for(uint i = 0; i < lambdaSize; i++) {
        outF[b] << std::setw(25) << fixStr << std::setprecision(4) <<
                  freeEnVal.lambdaCoulomb[i] << "," << std::setprecision(4) <<
                  freeEnVal.lambdaVDW[i] << ") ";
      }
      outF[b] << std::setw(25) << std::right << "PV(kJ/mol)\n";

      outF[b] << std::setprecision(std::numeric_limits<double>::digits10);
      outF[b].setf(std::ios_base::right, std::ios_base::adjustfield);
    } else
      std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                << "(Free Energy file)" << std::endl;
  }
}

void FreeEnergyOutput::CalculateFreeEnergy(const uint b)
{
  PV = var->pressure[b] * var->volumeRef[b] / var->numByBox[b];
  PV *= unit::K_TO_KJ_PER_MOL * unit::BAR_TO_K_MOLECULE_PER_A3;
  Etotal = var->energyRef[b].Total() / var->numByBox[b];
  Etotal *= unit::K_TO_KJ_PER_MOL;
  uint molIndex = lambdaRef.GetMolIndex(b);
  //Reset the energy value
  dUdL_VDW[b].Zero();
  dUdL_Coulomb[b].Zero();
  for(uint i = 0; i < lambdaSize; i++){
    energyDiff[b][i].Zero();
  }
  //Calculate delta E and dE/dLambda for all lambda states
  calcEn.EnergyChange(energyDiff[b], dUdL_VDW[b], dUdL_Coulomb[b], 
                      freeEnVal.lambdaVDW, freeEnVal.lambdaCoulomb, iState, 
                      molIndex, b);
  //Convert to kJ/mol
  dUdL_VDW[b] *= unit::K_TO_KJ_PER_MOL;
  dUdL_Coulomb[b] *= unit::K_TO_KJ_PER_MOL;
  for(uint i = 0; i < lambdaSize; i++){
    energyDiff[b][i] *= unit::K_TO_KJ_PER_MOL;
  }
}