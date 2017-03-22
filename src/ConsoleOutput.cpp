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
  if (step==0)
  {
    std::cout << std::endl << "########################## INITIAL SIMULATION ENERGY ##########################" << std::endl;
    for(uint b=0; b<BOX_TOTAL; b++)
    {
      PrintEnergyTitle(b);
      std::cout << std::endl;
      PrintEnergy(b, var->energyRef[b], var->virialRef[b]);
      std::cout << std::endl;
    }
    std::cout << "###############################################################################" << std::endl << std::endl;

    std::cout << "############################# STARTING SIMULATION #############################" << std::endl;
    for(uint b=0; b<BOX_TOTAL; b++)
    {
      PrintMoveTitle(b);
      std::cout << std::endl;
      if(enableEnergy)
      {
	PrintEnergyTitle(b);
      }
      PrintStatisticTitle(b);
      std::cout << std::endl;
    }
  }
  else
  {
    for(uint b=0; b<BOX_TOTAL; b++)
    {
      PrintMove(b, step);
      std::cout << std::endl;
      if(enableEnergy)
      {
	PrintEnergy(b, var->energyRef[b], var->virialRef[b]);
	std::cout << std::endl;
      }
      PrintStatistic(b);
      std::cout << std::endl;

      PrintPressureTensor(b);
      std::cout << std::endl;
    }
  }
}

void ConsoleOutput::PrintMove(const uint box, const ulong step) const
{
  uint sub;
  std::string title = "MOVE_BOX_";
  title += (box? "1": "0");
  printElement(title, elementWidth);

#if ENSEMBLE == GCMC
  if(box == mv::BOX0)
  {
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

#if ENSEMBLE == GEMC
  sub = mv::GetMoveSubIndex(mv::VOL_TRANSFER, box);
  printElement(var->GetTries(sub), elementWidth);
  printElement(var->GetAccepted(sub), elementWidth);
  printElement(var->GetAcceptPercent(sub), elementWidth);
  printElement(var->GetScale(sub), elementWidth);
#endif

  printElement(step + 1, elementWidth);
  std::cout << std::endl;
}

void ConsoleOutput::PrintStatistic(const uint box) const
{
  double density = 0.0;
  uint offset = box * var->numKinds;
  std::string title = "STAT_BOX_";
  title += (box? "1":"0");
  printElement(title, elementWidth);

  if(enableVolume)
    printElement(var->volumeRef[box] , elementWidth);
  if(enablePressure)
    printElement(var->pressure[box] , elementWidth);
  if(enableMol)
    printElement(var->numByBox[box], elementWidth);
  

  for(uint k=0; k<var->numKinds; k++)
  {
    uint kb = k+offset;
    if(var->numKinds > 1)
      printElement(var->molFractionByKindBox[kb], elementWidth);
  }

  if(enableDens)
    printElement(var->densityTot[box], elementWidth);
  if(enableSurfTension)
    printElement(var->surfaceTens[box], elementWidth);

  std::cout << std::endl;
}

void ConsoleOutput::PrintPressureTensor(const uint box) const
{
   std::string title = "PRESSURE_BOX_";
   title += (box? "1":"0");

   if(enablePressure)
   {
      printElement(title, elementWidth);
      for(uint i = 0; i < 3; i++)
      {
	 for(uint j = 0; j < 3; j++)
	 {
	   printElement(var->pressureTens[box][i][j], elementWidth);
	 }
      }
   }
   std::cout << std::endl;  
}


void ConsoleOutput::PrintEnergy(const uint box, Energy const& en, Virial const& vir) const
{
  std::string title = "ENERGY_BOX_";
  title += (box? "1":"0");
  printElement(title, elementWidth);

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

void ConsoleOutput::PrintEnergyTitle(const uint box)
{
  std::string title = "ETITLE_BOX_";
  title += (box? "1": "0");
  printElement(title, elementWidth);

  printElement("TOTAL", elementWidth);
  printElement("INTRA(B)", elementWidth);
  printElement("INTRA(N)", elementWidth);
  printElement("INTER(LJ)", elementWidth);
  printElement("LRC", elementWidth);
  printElement("TOTAL_ELECT", elementWidth);
  printElement("REAL", elementWidth);
  printElement("RECIP", elementWidth);
  printElement("SELF", elementWidth);
  printElement("CORR", elementWidth);
  std::cout << std::endl;
}

void ConsoleOutput::PrintStatisticTitle(const uint box)
{
  uint offset = box * var->numKinds;
  std::string title = "STITLE_BOX_";
  title += (box? "1": "0");
  printElement(title, elementWidth);

  if(enableVolume)
    printElement("VOLUME", elementWidth);
  if(enablePressure)
    printElement("PRESSURE", elementWidth);
  if(enableMol)
    printElement("TOTALMOL", elementWidth);
  

  for(uint k=0; k < var->numKinds; k++)
  {
    if(var->numKinds > 1)
    {
      std::string molName = "MOLFRAC_" + var->resKindNames[k];
      printElement(molName, elementWidth);
    }
  }

  if(enableDens)
    printElement("TOT_DENSITY", elementWidth);
  if(enableSurfTension)
    printElement("SURF_TENSION", elementWidth);

  std::cout << std::endl;
}


void ConsoleOutput::PrintMoveTitle(const uint box)
{
  std::string title = "MTITLE_BOX_";
  title += (box? "1": "0");
  printElement(title, elementWidth);

#if ENSEMBLE == GCMC
  if(box == mv::BOX0)
  {
#endif
    printElement("DISTRY", elementWidth);
    printElement("DISACCEPT", elementWidth);
    printElement("DISACCEPT%", elementWidth);
    printElement("DISMAX", elementWidth);

    printElement("ROTATE", elementWidth);
    printElement("ROTACCEPT", elementWidth);
    printElement("ROTACCEPT%", elementWidth);
    printElement("ROTMAX", elementWidth);

    printElement("INTRASWAP", elementWidth);
    printElement("INTACCEPT", elementWidth);
    printElement("INTACCEPT%", elementWidth);
    //printElement("INTMAX", elementWidth);
#if ENSEMBLE == GCMC
  }
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  printElement("TRANSFER", elementWidth);
  printElement("TRANACCEPT", elementWidth);
  printElement("TRANACCEPT%", elementWidth);
  //printElement("TRANMAX", elementWidth);
#endif

#if ENSEMBLE == GEMC
  printElement("VOLUME", elementWidth);
  printElement("VOLACCEPT", elementWidth);
  printElement("VOLACCEPT%", elementWidth);
  printElement("VOLMAX", elementWidth);
#endif

  printElement("STEP", elementWidth);
  std::cout << std::endl;
}

template <typename T>
void ConsoleOutput::printElement( const T t, const int width) const
{
  const char separator = ' ';
  std::cout << left << std::fixed << std::setprecision(4) << setw(width) <<
    setfill(separator) << t;
}
