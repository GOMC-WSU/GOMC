/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "ConsoleOutput.h" //For spec;

#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed

#include "EnsemblePreprocessor.h" //For BOX_TOTAL, ensemble
#include "FFConst.h"              //For density conv.
#include "MoleculeKind.h"         //For kind names
#include "MoveConst.h"            //For move index constants, name constants.
#include "OutputVars.h"
#include "PDBConst.h"   //For resname len.
#include "StaticVals.h" //for init
#include "System.h"     //for init

void ConsoleOutput::DoOutput(const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::CONSOLE_OUTPUT);
  if (WriteConsoleHeaders) {
    std::cout << std::endl
              << "#############################################################"
                 "###################"
              << std::endl;
    std::cout << "########################## INITIAL SIMULATION ENERGY "
                 "###########################"
              << std::endl
              << std::endl;

    PrintEnergyTitle();
    std::cout << std::endl;

    for (uint b = 0; b < BOX_TOTAL; b++) {
      PrintEnergy(b, var->energyRef[b], step - 1);
      std::cout << std::endl;
    }

    if (enableStat) {
      PrintStatisticTitle();
      std::cout << std::endl;

      for (uint b = 0; b < BOX_TOTAL; b++) {
        PrintStatistic(b, step - 1);
        std::cout << std::endl;
      }
    }

    std::cout << "#############################################################"
                 "###################"
              << std::endl;

    std::cout << "############################# STARTING SIMULATION "
                 "##############################"
              << std::endl
              << std::endl;

    if (!forceOutput) {
      PrintMoveTitle();
      std::cout << std::endl;
    }

    if (enableEnergy) {
      PrintEnergyTitle();
      std::cout << std::endl;
    }

    if (enableStat) {
      PrintStatisticTitle();
      std::cout << std::endl;
    }
    WriteConsoleHeaders = false;
  } else {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (!forceOutput) {
        PrintMove(b, step);
        std::cout << std::endl;
      }

      if (enableEnergy) {
        PrintEnergy(b, var->energyRef[b], step);
        std::cout << std::endl;
      }

      if (enableStat) {
        PrintStatistic(b, step);
        std::cout << std::endl;
      }

      if (enablePressure) {
        PrintPressureTensor(b, step);
        std::cout << std::endl;
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CONSOLE_OUTPUT);
}

void ConsoleOutput::DoOutputRestart(const ulong step) {}

void ConsoleOutput::PrintMove(const uint box, const ulong step) const {
  uint sub;
  std::string title = "MOVE_";
  title += (box ? "1:" : "0:");
  printElementStep(title, step + 1, elementWidth);

#if ENSEMBLE == GCMC
  if (box == mv::BOX0) {
#endif
    if (var->Performed(mv::DISPLACE)) {
      sub = mv::DISPLACE;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
      printElement(var->GetScale(box, sub), elementWidth);
    }

    if (var->Performed(mv::ROTATE)) {
      sub = mv::ROTATE;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
      printElement(var->GetScale(box, sub), elementWidth);
    }

    if (var->Performed(mv::MULTIPARTICLE)) {
      sub = mv::MULTIPARTICLE;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::MULTIPARTICLE_BM)) {
      sub = mv::MULTIPARTICLE_BM;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::INTRA_SWAP)) {
      sub = mv::INTRA_SWAP;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::REGROWTH)) {
      sub = mv::REGROWTH;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::INTRA_MEMC)) {
      sub = mv::INTRA_MEMC;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::CRANKSHAFT)) {
      sub = mv::CRANKSHAFT;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

    if (var->Performed(mv::INTRA_TARGETED_SWAP)) {
      sub = mv::INTRA_TARGETED_SWAP;
      printElement(var->GetTries(box, sub), elementWidth);
      printElement(var->GetAccepted(box, sub), elementWidth);
      printElement(var->GetAcceptPercent(box, sub), elementWidth);
    }

#if ENSEMBLE == GCMC
  }
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  if (var->Performed(mv::MOL_TRANSFER)) {
    sub = mv::MOL_TRANSFER;
    printElement(var->GetTries(box, sub), elementWidth);
    printElement(var->GetAccepted(box, sub), elementWidth);
    printElement(var->GetAcceptPercent(box, sub), elementWidth);
  }

  if (var->Performed(mv::MEMC)) {
    sub = mv::MEMC;
    printElement(var->GetTries(box, sub), elementWidth);
    printElement(var->GetAccepted(box, sub), elementWidth);
    printElement(var->GetAcceptPercent(box, sub), elementWidth);
  }

  if (var->Performed(mv::NE_MTMC)) {
    sub = mv::NE_MTMC;
    printElement(var->GetTries(box, sub), elementWidth);
    printElement(var->GetAccepted(box, sub), elementWidth);
    printElement(var->GetAcceptPercent(box, sub), elementWidth);
  }

  if (var->Performed(mv::TARGETED_SWAP)) {
    sub = mv::TARGETED_SWAP;
    printElement(var->GetTries(box, sub), elementWidth);
    printElement(var->GetAccepted(box, sub), elementWidth);
    printElement(var->GetAcceptPercent(box, sub), elementWidth);
  }
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  if (var->Performed(mv::VOL_TRANSFER)) {
    sub = mv::VOL_TRANSFER;
    printElement(var->GetTries(box, sub), elementWidth);
    printElement(var->GetAccepted(box, sub), elementWidth);
    printElement(var->GetAcceptPercent(box, sub), elementWidth);
    printElement(var->GetScale(box, sub), elementWidth);
  }
#endif

  std::cout << std::endl;
}

void ConsoleOutput::PrintStatistic(const uint box, const ulong step) const {
  uint offset = box * var->numKinds;

  std::string title = "STAT_";
  sstrm::Converter toStr;
  std::string numStr = "";
  toStr << box;
  toStr >> numStr;
  title += numStr + ":";
  printElementStep(title, step + 1, elementWidth);

  if (enableVolume)
    printElement(var->volumeRef[box], elementWidth);

  if (enablePressure) {
    printElement(var->pressure[box], elementWidth);
    printElement(var->compressibility[box], elementWidth);
    printElement(var->enthalpy[box], elementWidth);
#if ENSEMBLE == GEMC
    printElement(var->heatOfVap, elementWidth);
#endif
  }
  if (enableMol) {
    printElement(var->numByBox[box], elementWidth);

    for (uint k = 0; k < var->numKinds; k++) {
      uint kb = k + offset;
      if (var->numKinds > 1)
        printElement(var->molFractionByKindBox[kb], elementWidth, 8);
    }
    for (uint k = 0; k < var->numKinds; k++) {
      uint kb = k + offset;
      if (var->numKinds > 1)
        printElement(var->densityByKindBox[kb], elementWidth, 8);
    }
  }

  if (enableDens)
    printElement(var->densityTot[box], elementWidth);
  if (enableSurfTension)
    printElement(var->surfaceTens[box], elementWidth);

  std::cout << std::endl;
}

void ConsoleOutput::PrintPressureTensor(const uint box,
                                        const ulong step) const {
  std::string title = "PRES_";
  sstrm::Converter toStr;
  std::string numStr = "";
  toStr << box;
  toStr >> numStr;
  title += numStr + ":";
  printElementStep(title, step + 1, elementWidth);

  for (uint i = 0; i < 3; i++) {
    // If you calculate the pressure tensor for W12, W13, W23 we print all 9
    // values of tensor
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

void ConsoleOutput::PrintEnergy(const uint box, Energy const &en,
                                const ulong step) const {
  std::string title = "ENER_";
  sstrm::Converter toStr;
  std::string numStr = "";
  toStr << box;
  toStr >> numStr;
  title += numStr + ":";
  printElementStep(title, step + 1, elementWidth);

  printElement(en.total, elementWidth);
  printElement(en.intraBond, elementWidth);
  printElement(en.intraNonbond, elementWidth);
  printElement(en.inter, elementWidth);
  printElement(en.tailCorrection, elementWidth);
  printElement(en.totalElect, elementWidth);
  printElement(en.real, elementWidth);
  printElement(en.recip, elementWidth);
  printElement(en.self, elementWidth);
  printElement(en.correction, elementWidth);
  std::cout << std::endl;
}

void ConsoleOutput::PrintEnergyTitle() {
  std::string title = "ETITLE:";
  size_t loc = title.length();
  title += "STEP";
  title.insert(loc, elementWidth - title.length(), ' ');
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

void ConsoleOutput::PrintStatisticTitle() {
  // uint offset = box * var->numKinds;
  if (enableStat) {
    std::string title = "STITLE:";
    size_t loc = title.length();
    title += "STEP";
    title.insert(loc, elementWidth - title.length(), ' ');
    printElement(title, elementWidth);
  }

  if (enableVolume)
    printElement("VOLUME", elementWidth);

  if (enablePressure) {
    printElement("PRESSURE", elementWidth);
    printElement("COMPRESSIBILITY", elementWidth);
    printElement("ENTHALPY", elementWidth);
#if ENSEMBLE == GEMC
    printElement("HEAT_VAP", elementWidth);
#endif
  }

  if (enableMol) {
    printElement("TOTALMOL", elementWidth);

    for (uint k = 0; k < var->numKinds; k++) {
      if (var->numKinds > 1) {
        std::string molName = "MOLFRAC_" + var->molKindNames[k];
        printElement(molName, elementWidth);
      }
    }
    for (uint k = 0; k < var->numKinds; k++) {
      if (var->numKinds > 1) {
        std::string molName = "MOLDENS_" + var->molKindNames[k];
        printElement(molName, elementWidth);
      }
    }
  }

  if (enableDens)
    printElement("TOT_DENSITY", elementWidth);
  if (enableSurfTension)
    printElement("SURF_TENSION", elementWidth);

  std::cout << std::endl;
}

void ConsoleOutput::PrintMoveTitle() {
  std::string title = "MTITLE:";
  size_t loc = title.length();
  title += "STEP";
  title.insert(loc, elementWidth - title.length(), ' ');
  printElement(title, elementWidth);
  if (var->Performed(mv::DISPLACE)) {
    printElement("DISTRY", elementWidth);
    printElement("DISACCEPT", elementWidth);
    printElement("DISACCEPT%", elementWidth);
    printElement("DISMAX", elementWidth);
  }

  if (var->Performed(mv::ROTATE)) {
    printElement("ROTATE", elementWidth);
    printElement("ROTACCEPT", elementWidth);
    printElement("ROTACCEPT%", elementWidth);
    printElement("ROTMAX", elementWidth);
  }

  if (var->Performed(mv::MULTIPARTICLE)) {
    printElement("MULTIPARTICLE", elementWidth);
    printElement("MPACCEPT", elementWidth);
    printElement("MPACCEPT%", elementWidth);
  }

  if (var->Performed(mv::MULTIPARTICLE_BM)) {
    printElement("MULTIPARTICLEBM", elementWidth);
    printElement("MPBMACCEPT", elementWidth);
    printElement("MPBMACCEPT%", elementWidth);
  }

  if (var->Performed(mv::INTRA_SWAP)) {
    printElement("INTRASWAP", elementWidth);
    printElement("INTACCEPT", elementWidth);
    printElement("INTACCEPT%", elementWidth);
  }

  if (var->Performed(mv::REGROWTH)) {
    printElement("REGROWTH", elementWidth);
    printElement("REGROWACCEPT", elementWidth);
    printElement("REGROWACCEPT%", elementWidth);
  }

  if (var->Performed(mv::INTRA_MEMC)) {
    printElement("INTRAMOLEXCHANGE", elementWidth);
    printElement("INTMOLEXCACCEPT", elementWidth);
    printElement("INTMOLEXACCEPT%", elementWidth);
  }

  if (var->Performed(mv::CRANKSHAFT)) {
    printElement("CRANKSHAFT", elementWidth);
    printElement("CRKSHAFTACCEPT", elementWidth);
    printElement("CRKSHAFTACCEPT%", elementWidth);
  }

  if (var->Performed(mv::INTRA_TARGETED_SWAP)) {
    printElement("INTRATARGET", elementWidth);
    printElement("INTARGACCEPT", elementWidth);
    printElement("INTARGACCEPT%", elementWidth);
  }

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  if (var->Performed(mv::MOL_TRANSFER)) {
    printElement("TRANSFER", elementWidth);
    printElement("TRANACCEPT", elementWidth);
    printElement("TRANACCEPT%", elementWidth);
  }

  if (var->Performed(mv::MEMC)) {
    printElement("MOLEXCHANGE", elementWidth);
    printElement("MOLEXACCEPT", elementWidth);
    printElement("MOLEXACCEPT%", elementWidth);
  }

  if (var->Performed(mv::NE_MTMC)) {
    printElement("NEMTMCTRANSF", elementWidth);
    printElement("NEMTMCACCEPT", elementWidth);
    printElement("NEMTMCACCEPT%", elementWidth);
  }

  if (var->Performed(mv::TARGETED_SWAP)) {
    printElement("TARGETTRANSFER", elementWidth);
    printElement("TARGETTACCEPT", elementWidth);
    printElement("TARGETTACCEPT%", elementWidth);
  }
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  if (var->Performed(mv::VOL_TRANSFER)) {
    printElement("VOLUME", elementWidth);
    printElement("VOLACCEPT", elementWidth);
    printElement("VOLACCEPT%", elementWidth);
    printElement("VOLMAX", elementWidth);
  }
#endif

  std::cout << std::endl;
}

void ConsoleOutput::printElement(const double t, const int width,
                                 uint precision) const {
  const char separator = ' ';
  if (std::abs(t) > 1e99) {
    std::cout << std::right << std::scientific
              << std::setprecision(precision - 1) << std::setw(width)
              << std::setfill(separator) << t;
  } else {
    std::cout << std::right << std::scientific << std::setprecision(precision)
              << std::setw(width) << std::setfill(separator) << t;
  }
}

void ConsoleOutput::printElement(const uint t, const int width) const {
  const char separator = ' ';
  std::cout << std::right << std::scientific << std::setw(width)
            << std::setfill(separator) << t;
}

void ConsoleOutput::printElement(const std::string t, const int width) const {
  const char separator = ' ';
  std::cout << std::right << std::scientific << std::setw(width)
            << std::setfill(separator) << t;
}

template <typename T>
void ConsoleOutput::printElementStep(const T t, const ulong step,
                                     const int width) const {
  std::cout << t << std::right << std::setw(width - 7) << step;
}
