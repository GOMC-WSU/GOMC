/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "ConsoleOutput.h"          //For spec;
#include "EnsemblePreprocessor.h"   //For BOX_TOTAL, ensemble
#include "MoveConst.h"    //For move index constants, name constants.
#include "FFConst.h"                //For density conv.
#include "System.h"                 //for init
#include "StaticVals.h"             //for init
#include "MoleculeKind.h"           //For kind names
#include "PDBConst.h"               //For resname len.
#include "OutputVars.h"

#include <iostream>                 // std::cout, std::fixed
#include <iomanip>                  // std::setprecision

void ConsoleOutput::DoOutput(const ulong step)
{
  if (step == 0) {
    std::cout << std::endl << "################################################################################" << std::endl;
    std::cout << "########################## INITIAL SIMULATION ENERGY ###########################" << std::endl << std::endl;

    PrintEnergyTitle();
    std::cout << std::endl;

    for (uint b = 0; b < BOX_TOTAL; b++) {
      PrintEnergy(b, var->energyRef[b], var->virialRef[b], -1);
      std::cout <<  std::endl;
    }

    if(enableStat) {
      PrintStatisticTitle();
      std::cout << std::endl;

      for (uint b = 0; b < BOX_TOTAL; b++) {
        PrintStatistic(b, -1);
        std::cout << std::endl;
      }
    }

    std::cout << "################################################################################" << std::endl;

    std::cout << "############################# STARTING SIMULATION ##############################" << std::endl << std::endl;

    PrintMoveTitle();
    std::cout << std::endl;

    if(enableEnergy) {
      PrintEnergyTitle();
      std::cout << std::endl;
    }

    if(enableStat) {
      PrintStatisticTitle();
      std::cout << std::endl;
    }
  } else {
    for(uint b = 0; b < BOX_TOTAL; b++) {
      PrintMove(b, step);
      std::cout << std::endl;

      if(enableEnergy) {
        PrintEnergy(b, var->energyRef[b], var->virialRef[b], step);
        std::cout << std::endl;
      }

      if(enableStat) {
        PrintStatistic(b, step);
        std::cout << std::endl;
      }

      if(enablePressure) {
        PrintPressureTensor(b, step);
        std::cout << std::endl;
      }

    }

  }
}

void ConsoleOutput::PrintMove(const uint box, const ulong step) const
{
  uint sub;
  std::string title = "MOVE_";
  title += (box ? "1:" : "0:");
  printElementStep(title, step + 1, elementWidth);

#if ENSEMBLE == GCMC
  if(box == mv::BOX0) {
#endif
    sub = mv::GetMoveSubIndex(mv::DISPLACE, box);
    printElement(var->GetTries(sub), elementWidth);
    printElement(var->GetAccepted(sub), elementWidth);
    printElement(var->GetAcceptPercent(sub), elementWidth);
    printElement(var->GetScale(sub), elementWidth);

    sub = mv::GetMoveSubIndex(mv::ROTATE, box);
    printElement(var->GetTries(sub), elementWidth);
    printElement(var->GetAccepted(sub), elementWidth);
    printElement(var->GetAcceptPercent(sub), elementWidth);
    printElement(var->GetScale(sub), elementWidth);

    sub = mv::GetMoveSubIndex(mv::MULTIPARTICLE, box);
    printElement(var->GetTries(sub), elementWidth);
    printElement(var->GetAccepted(sub), elementWidth);
    printElement(var->GetAcceptPercent(sub), elementWidth);
    printElement(var->GetScale(sub), elementWidth);

    sub = mv::GetMoveSubIndex(mv::INTRA_SWAP, box);
    printElement(var->GetTries(sub), elementWidth);
    printElement(var->GetAccepted(sub), elementWidth);
    printElement(var->GetAcceptPercent(sub), elementWidth);
    //printElement(var->GetScale(sub), elementWidth);
#if ENSEMBLE == GCMC
  }
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  sub = mv::GetMoveSubIndex(mv::MOL_TRANSFER, box);
  printElement(var->GetTries(sub), elementWidth);
  printElement(var->GetAccepted(sub), elementWidth);
  printElement(var->GetAcceptPercent(sub), elementWidth);
  //printElement(var->GetScale(sub), elementWidth);
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  sub = mv::GetMoveSubIndex(mv::VOL_TRANSFER, box);
  printElement(var->GetTries(sub), elementWidth);
  printElement(var->GetAccepted(sub), elementWidth);
  printElement(var->GetAcceptPercent(sub), elementWidth);
  printElement(var->GetScale(sub), elementWidth);
#endif

  std::cout << std::endl;
}

void ConsoleOutput::PrintStatistic(const uint box, const ulong step) const
{
  double density = 0.0;
  uint offset = box * var->numKinds;

  std::string title = "STAT_";
  title += (box ? "1:" : "0:");
  printElementStep(title, step + 1, elementWidth);

  if(enableVolume)
    printElement(var->volumeRef[box], elementWidth);

  if(enablePressure)
    printElement(var->pressure[box], elementWidth);

  if(enableMol) {
    printElement(var->numByBox[box], elementWidth);

    for(uint k = 0; k < var->numKinds; k++) {
      uint kb = k + offset;

      if(var->numKinds > 1)
        printElement(var->molFractionByKindBox[kb], elementWidth, 6);
    }
  }

  if(enableDens)
    printElement(var->densityTot[box], elementWidth);
  if(enableSurfTension)
    printElement(var->surfaceTens[box], elementWidth);

  std::cout << std::endl;
}

void ConsoleOutput::PrintPressureTensor(const uint box, const ulong step) const
{
  std::string title = "PRES_";
  title += (box ? "1:" : "0:");
  printElementStep(title, step + 1, elementWidth);

  for(uint i = 0; i < 3; i++) {
    //If you calculate the pressure tensor for W12, W13, W23 we print all 9 values of tensor
    /*
      for(uint j = 0; j < 3; j++)
      {
         printElement(var->pressureTens[box][i][j], elementWidth);
      }
    */
    // Else, just print the diameter of the pressure Tensor, W11, W22, W33
    printElement(var->pressureTens[box][i][i], elementWidth);
  }
  std::cout << std::endl;
}


void ConsoleOutput::PrintEnergy(const uint box, Energy const& en,
                                Virial const& vir, const ulong step) const
{
  std::string title = "ENER_";
  title += (box ? "1:" : "0:");
  printElementStep(title, step + 1, elementWidth);

  printElement(en.total, elementWidth);
  printElement(en.intraBond, elementWidth);
  printElement(en.intraNonbond, elementWidth);
  printElement(en.inter, elementWidth);
  printElement(en.tc, elementWidth);
  printElement(en.totalElect, elementWidth);
  printElement(en.real, elementWidth);
  printElement(en.recip, elementWidth);
  printElement(en.self, elementWidth);
  printElement(en.correction, elementWidth);
  std::cout << std::endl;
}

void ConsoleOutput::PrintEnergyTitle()
{
  std::string title = "ETITLE:";
  title += "     STEP";
  printElement(title, elementWidth);

  printElement("TOTAL", elementWidth);
  printElement("INTRA(B)", elementWidth);
  printElement("INTRA(NB)", elementWidth);
  printElement("INTER(LJ)", elementWidth);
  printElement("LRC", elementWidth);
  printElement("TOTAL_ELECT", elementWidth);
  printElement("REAL", elementWidth);
  printElement("RECIP", elementWidth);
  printElement("SELF", elementWidth);
  printElement("CORR", elementWidth);
  std::cout << std::endl;
}

void ConsoleOutput::PrintStatisticTitle()
{
  //uint offset = box * var->numKinds;
  if(enableStat) {
    std::string title = "STITLE:";
    title += "     STEP";
    printElement(title, elementWidth);
  }

  if(enableVolume)
    printElement("VOLUME", elementWidth);

  if(enablePressure)
    printElement("PRESSURE", elementWidth);

  if(enableMol) {
    printElement("TOTALMOL", elementWidth);

    for(uint k = 0; k < var->numKinds; k++) {
      if(var->numKinds > 1) {
        std::string molName = "MOLFRAC_" + var->resKindNames[k];
        printElement(molName, elementWidth);
      }
    }
  }

  if(enableDens)
    printElement("TOT_DENSITY", elementWidth);
  if(enableSurfTension)
    printElement("SURF_TENSION", elementWidth);

  std::cout << std::endl;
}


void ConsoleOutput::PrintMoveTitle()
{
  std::string title = "MTITLE:";
  title += "     STEP";
  printElement(title, elementWidth);

  printElement("DISTRY", elementWidth);
  printElement("DISACCEPT", elementWidth);
  printElement("DISACCEPT%", elementWidth);
  printElement("DISMAX", elementWidth);

  printElement("ROTATE", elementWidth);
  printElement("ROTACCEPT", elementWidth);
  printElement("ROTACCEPT%", elementWidth);
  printElement("ROTMAX", elementWidth);

  printElement("MULTIPARTICLE", elementWidth);
  printElement("MPACCEPT", elementWidth);
  printElement("MPACCEPT%", elementWidth);
  printElement("MPMAX", elementWidth);

  printElement("INTRASWAP", elementWidth);
  printElement("INTACCEPT", elementWidth);
  printElement("INTACCEPT%", elementWidth);


#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  printElement("TRANSFER", elementWidth);
  printElement("TRANACCEPT", elementWidth);
  printElement("TRANACCEPT%", elementWidth);
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  printElement("VOLUME", elementWidth);
  printElement("VOLACCEPT", elementWidth);
  printElement("VOLACCEPT%", elementWidth);
  printElement("VOLMAX", elementWidth);
#endif

  std::cout << std::endl;
}

void ConsoleOutput::printElement(const double t, const int width,
                                 uint percision) const
{
  const char separator = ' ';
  if(abs(t) > 9999999999.9999) {
    std::cout << right << std::fixed << std::setprecision(0) <<
              setw(width) << setfill(separator) << 9999999999;
  } else {
    std::cout << right << std::fixed << std::setprecision(percision) <<
              setw(width) << setfill(separator) << t;
  }

}

void ConsoleOutput::printElement(const uint t, const int width) const
{
  const char separator = ' ';
  std::cout << right << std::fixed  << setw(width) <<
            setfill(separator) << t;
}

void ConsoleOutput::printElement(const std::string t, const int width) const
{
  const char separator = ' ';
  std::cout << right << std::fixed << setw(width) <<
            setfill(separator) << t;
}

template <typename T>
void ConsoleOutput::printElementStep( const T t, const ulong step,
                                      const int width) const
{
  std::cout << t << right << setw(width - 7) << step;
}
