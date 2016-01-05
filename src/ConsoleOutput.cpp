#include "ConsoleOutput.h"          //For spec;
#include "EnsemblePreprocessor.h"   //For BOX_TOTAL, ensemble
#include "MoveConst.h"    //For move index constants, name constants.
#include "FFConst.h"                //For density conv.
#include "System.h"                 //for init
#include "StaticVals.h"             //for init  
#include "MoleculeKind.h"           //For kind names
#include "PDBConst.h"               //For resname len.
#include "OutputVars.h"

#include <iostream>                 // for cout;

void ConsoleOutput::DoOutput(const ulong step)
{
   if (step==0)
      std::cout << "STARTING SIMULATION!!" << std::endl << std::endl;
   else 
      std::cout << "STEP " << step + 1 << " of " << totSimSteps << " ("
		<< (double)(step + 1)/(double)(totSimSteps) * 100.0
		<< "% done)" << std::endl << std::endl;

   for (uint b = 0; b < BOX_TOTAL; b++)
      PrintBox(b, step);

   PrintSysStat();
}

void ConsoleOutput::PrintBox(const uint box, const ulong step) const
{
   std::cout << std::endl;
   std::string banner = "BOX ";
   uint offset = box * var->numKinds;
   banner += (box?"1 ":"0 ");
   PrintBanner(banner);
   
   std::cout << "Volume: " << var->volumeRef[box] << " A^3" << std::endl
	     << "Tot # mol: " << var->numByBox[box] << std::endl << std::endl;
   
   for(uint k = 0; k < var->numKinds; ++k)
      PrintMolKind(k, k+offset);

   PrintEnergy(var->energyRef[box], var->virialRef[box],
               (box < BOXES_WITH_U_NB));

   if (step > 0 )
   {
      PrintMoveStat(box, step);
   }
}

void ConsoleOutput::PrintSysStat() const
{
   PrintBanner("SYSTEM");
   PrintEnergy(*(var->energyTotRef), *(var->virialTotRef));
}

void ConsoleOutput::PrintMoveKind(bool & somethingPrinted,
				  const uint m,
				  const uint b,
				  const ulong step) const
{      
   uint sub = mv::GetMoveSubIndex(m, b);
   if (var->movePercRef[m] == 0.0 && step==0)
   {
      somethingPrinted |= true;
      std::cout << "Reminder: " << mv::MOVE_NAME[sub] 
		<< " move is off!" << std::endl;
   }
   else if (var->movePercRef[m] != 0.0 && step != 0)
   {
      somethingPrinted |= true;
      std::cout << mv::MOVE_NAME[sub] << " -- tries: " << var->GetTries(sub)
		<< "; # Accept: " << var->GetAccepted(sub)
		<< "; % Accept : " << var->GetAcceptPercent(sub);
      if (sub < mv::SCALEABLE)
	 std::cout << "; Max Amt.: " << var->GetScale(sub);
      std::cout << std::endl;
   }
}

void ConsoleOutput::PrintMoveStat(const uint box, const ulong step) const
{
   bool somethingPrinted = false;
#if ENSEMBLE == GCMC
   if (box == mv::BOX0)
   {
#endif
      PrintMoveKind(somethingPrinted, mv::DISPLACE, box, step);
      PrintMoveKind(somethingPrinted, mv::ROTATE, box, step);
      PrintMoveKind(somethingPrinted, mv::INTRA_SWAP, box, step);
#if ENSEMBLE == GCMC
   }
#endif

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
   PrintMoveKind(somethingPrinted, mv::MOL_TRANSFER, box, step);
#endif

#if ENSEMBLE == GEMC
   PrintMoveKind(somethingPrinted, mv::VOL_TRANSFER, box, step);
#endif
   if (somethingPrinted)
      std::cout << std::endl;
}

void ConsoleOutput::PrintMolKind(const uint k, const uint kb) const
{
   std::cout << "Molecule Type " << var->resKindNames[k] << std::endl
	     << "-------------------" << std::endl
	     << "Number: " << var->numByKindBox[kb] << " ;" << std::endl;
   if (k > 1)
      std::cout <<  "Mol. Fraction: " << var->molFractionByKindBox[kb]
                << " ;" << std::endl;
   std::cout << "Density: "  << var->densityByKindBox[kb] 
	     << " g/ml" << std::endl << std::endl;
}

void ConsoleOutput::PrintEnergy(Energy const& en, Virial const& vir,
                                const bool intraOnly) const
{
   if (intraOnly)
      std::cout << "Energy (in K): total: " << en.total << std::endl
		<< "intra (bonded): " << en.intraBond << std::endl
		<< "intra (nonbonded): " << en.intraNonbond << std::endl
	        << "inter: " << en.inter << std::endl
		<< "inter (tail corr.): " << en.tc<< std::endl
                << "totalElect: " << en.totalElect << std::endl
		<< "real: " << en.real << std::endl
		<< "recip: " << en.recip << std::endl
		<< "self: " << en.self << std::endl
		<< "correction: " << en.correction << std::endl;
   /*      std::cout << "Energy (in K): total: " << en.total << std::endl
                << "inter: " << en.inter << " ; inter (tail corr.): " 
                << en.tc << std::endl << "intra (bonded): "
                << en.intraBond << " ; intra (nonbonded): " << en.intraNonbond 
				 << std::endl
		 << "; elect: " << en.totalElect  <<"; real: " << en.real
		<< "; recip: " << en.recip  << "; self: " << en.self
		<< "; correction: " << en.correction  << std::endl;
   */
}

void ConsoleOutput::PrintBanner(std::string const& str) const 
{
   std::cout << "------------------------" << std::endl
	     << "--    ==========      --" << std::endl
	     << "--  === "<< str << " ===    --" << std::endl
	     << "--    ==========      --" << std::endl
	     << "------------------------" << std::endl << std::endl;
}

