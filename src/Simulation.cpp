#include "Simulation.h"
#include "Setup.h"          //For setup object

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>

Simulation::Simulation(char const*const configFileName)
{
   //NOTE:
   //IMPORTANT! Keep this order...
   //as system depends on staticValues, and cpu sometimes depends on both.
   Setup set;
   set.Init(configFileName);
   totalSteps = set.config.sys.step.total;
   staticValues = new StaticVals();
   system = new System(*staticValues);
   staticValues->Init(set, *system); 
   system->Init(set);
   cpu = new CPUSide(*system, *staticValues);
   cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
             totalSteps);

   //Dump combined PSF
   PSFOutput psfOut(staticValues->mol, set.mol.kindMap, 
		    set.pdb.atoms.resKindNames);
   psfOut.PrintPSF(set.config.out.state.files.psf.name);
   std::cout << "Printed combined psf to file " 
	     << set.config.out.state.files.psf.name << '\n';

}

Simulation::~Simulation()
{
   delete staticValues;
   delete system;
   delete cpu;
}



void Simulation::RunSimulation(void)
{
   for (ulong step = 0; step < totalSteps; step++)
   {
      system->ChooseAndRunMove(step);
      cpu->Output(step);
#ifndef NDEBUG
      if ((step + 1) % 100 == 0)
         RunningCheck(step);
#endif
   }
}

#ifndef NDEBUG
void Simulation::RunningCheck(const uint step)
{
   SystemPotential pot = system->calcEnergy.SystemTotal();
   
   std::cout 
      << "================================================================="
      << std::endl << "-------------------------" << std::endl
      << " STEP: " << step
      << std::endl << "-------------------------" << std::endl
      << "Energy       INTRA B |     INTRA NB |        INTER |           TC |         REAL |         SELF |   CORRECTION |        RECIP"
      << std::endl
      << "System: "
      << std::setw(12) << system->potential.totalEnergy.intraBond << " | "
      << std::setw(12) << system->potential.totalEnergy.intraNonbond << " | "
      << std::setw(12) << system->potential.totalEnergy.inter << " | "
      << std::setw(12) << system->potential.totalEnergy.tc << " | "
      << std::setw(12) << system->potential.totalEnergy.real << " | "
      << std::setw(12) << system->potential.totalEnergy.self << " | "
      << std::setw(12) << system->potential.totalEnergy.correction << " | "
      << std::setw(12) << system->potential.totalEnergy.recip << std::endl
      << "Recalc: "
      << std::setw(12) << pot.totalEnergy.intraBond << " | "
      << std::setw(12) << pot.totalEnergy.intraNonbond << " | "
      << std::setw(12) << pot.totalEnergy.inter << " | "
      << std::setw(12) << pot.totalEnergy.tc << " | "
      << std::setw(12) << pot.totalEnergy.real << " | "
      << std::setw(12) << pot.totalEnergy.self << " | "
      << std::setw(12) << pot.totalEnergy.correction << " | "
      << std::setw(12) << pot.totalEnergy.recip << std::endl
      << "-------------------------" << std::endl
      << "Virial            INTER |           TC" << std::endl
      << "System: "
      << std::setw(15) << system->potential.totalVirial.inter << " | "
      << std::setw(12) << system->potential.totalVirial.tc << std::endl
      << "Recalc: "
      << std::setw(15) << pot.totalVirial.inter << " | "
      << std::setw(12) << pot.totalVirial.tc << std::endl
      << "-------------------------" << std::endl
      << "================================================================"
      << std::endl << std::endl;

}
#endif
