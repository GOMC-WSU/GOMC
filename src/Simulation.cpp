/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object

#include "EnergyTypes.h"
#include "PSFOutput.h"
#include <iostream>
#include <iomanip>

Simulation::Simulation(char const*const configFileName, MultiSim const*const& multisim):ms(multisim)
{
  startStep = 0;
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  set.Init(configFileName, multisim);
  totalSteps = set.config.sys.step.total;
  staticValues = new StaticVals(set);
  system = new System(*staticValues);
  staticValues->Init(set, *system);
  system->Init(set, startStep);
  //recal Init for static value for initializing ewald since ewald is
  //initialized in system
  staticValues->InitOver(set, *system);
  cpu = new CPUSide(*system, *staticValues);
  cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
            totalSteps, startStep);

  //Dump combined PSF
  PSFOutput psfOut(staticValues->mol, *system, set.mol.kindMap,
                   set.pdb.atoms.resKindNames);
  psfOut.PrintPSF(set.config.out.state.files.psf.name);
  std::cout << "Printed combined psf to file "
            << set.config.out.state.files.psf.name << '\n';

  if(totalSteps == 0) {
    frameSteps = set.pdb.GetFrameSteps(set.config.in.files.pdb.name);
  }
#if GOMC_LIB_MPI
  PTUtils = set.config.sys.step.parallelTemp ? new ParallelTemperingUtilities(ms, *system, *staticValues, set.config.sys.step.parallelTempFreq): NULL;
  exchangeResults.resize(ms->worldSize, false);
#endif
}

Simulation::~Simulation()
{
  delete cpu;
  delete system;
  delete staticValues;
}

void Simulation::RunSimulation(void)
{
  double startEnergy = system->potential.totalEnergy.total;
  if(totalSteps == 0) {
    for(int i = 0; i < frameSteps.size(); i++) {
      if(i == 0) {
        cpu->Output(frameSteps[0] - 1);
        continue;
      }
      system->RecalculateTrajectory(set, i + 1);
      cpu->Output(frameSteps[i] - 1);
    }
  }
  for (ulong step = startStep; step < totalSteps; step++) {
    system->moveSettings.AdjustMoves(step);
    system->ChooseAndRunMove(step);
    cpu->Output(step);

    if((step + 1) == cpu->equilSteps) {
      double currEnergy = system->potential.totalEnergy.total;
      if(abs(currEnergy - startEnergy) > 1.0e+10) {
        printf("Info: Recalculating the total energies to insure the accuracy"
               " of the computed \n"
               "      running energies.\n\n");
        system->calcEwald->Init();
        system->potential = system->calcEnergy.SystemTotal();
      }
    }
  #if GOMC_LIB_MPI

    if( staticValues->simEventFreq.parallelTemp && step > cpu->equilSteps && step % staticValues->simEventFreq.parallelTempFreq == 0){
      system->potential = system->calcEnergy.SystemTotal();
      exchangeResults = PTUtils->evaluateExchangeCriteria(step);
    
      //if (ms->worldRank == 1){
      if (exchangeResults[ms->worldRank] == true){


        std::cout << "A swap took place" << std::endl;

        //SystemPotential myPotentialCloneBeforeExchange(system->potential);
        //SystemPotential potBuffer(system->potential);

        PTUtils->exchangePositions(system->coordinates, ms, ms->worldRank-1, true);
        PTUtils->exchangeCOMs(system->com, ms, ms->worldRank-1, true);

     //   PTUtils->exchangeCellLists(myCellListCloneBeforeExchange, ms, ms->worldRank-1, true);
        PTUtils->exchangeCellLists(system->cellList, ms, ms->worldRank-1, true);
        //PTUtils->exchangePotentials(myPotentialCloneBeforeExchange, ms, ms->worldRank-1, true);
        PTUtils->exchangePotentials(system->potential, ms, ms->worldRank-1, true);

     //   system->cellList.GridAll(system->boxDimRef, system->coordinates, system->molLookup);
 //       system->potential = system->calcEnergy.SystemTotal();
        //potBuffer = system->calcEnergy.SystemTotal();
/*
        if (!potBuffer.ComparePotentials(myPotentialCloneBeforeExchange)){
          std::cout << "Potential objects have different states. Exiting!" << std::endl;
          exit(EXIT_FAILURE);
        } else {
          std::cout << "Potential objects have equivalent states." << std::endl;
        }
*/
      


/*
        if (!system->cellList.CompareCellList(myCellListCloneBeforeExchange, system->coordinates.Count())){
          std::cout << "Cell List objects have different states, not simply different orders. Exiting!" << std::endl;
          exit(EXIT_FAILURE);
        } else {
          std::cout << "Cell List objects have equivalent states." << std::endl;
        }
*/
       // system->calcEwald->UpdateVectorsAndRecipTerms();


//      } else if (ms->worldRank == 0){
      } else if(ms->worldRank+1 != ms->worldSize && exchangeResults[ms->worldRank+1] == true) {

        std::cout << "A swap took place" << std::endl;

        //SystemPotential myPotentialCloneBeforeExchange(system->potential);
        //SystemPotential potBuffer(system->potential);

        PTUtils->exchangePositions(system->coordinates, ms, ms->worldRank+1, false);
        PTUtils->exchangeCOMs(system->com, ms, ms->worldRank+1, false);

        //PTUtils->exchangeCellLists(myCellListCloneBeforeExchange, ms, ms->worldRank-1, true);
        PTUtils->exchangeCellLists(system->cellList, ms, ms->worldRank+1, false);
        PTUtils->exchangePotentials(system->potential, ms, ms->worldRank+1, false);

     //   system->cellList.GridAll(system->boxDimRef, system->coordinates, system->molLookup);
        //system->potential = system->calcEnergy.SystemTotal();
        //potBuffer = system->calcEnergy.SystemTotal();
/*
        if (!potBuffer.ComparePotentials(myPotentialCloneBeforeExchange)){
          std::cout << "Potential objects have different states. Exiting!" << std::endl;
          exit(EXIT_FAILURE);
        } else {
          std::cout << "Potential objects have equivalent states." << std::endl;
        }
*/
      }

    }
  #endif

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  system->PrintAcceptance();
  system->PrintTime();
}

#ifndef NDEBUG
void Simulation::RunningCheck(const uint step)
{
  system->calcEwald->UpdateVectorsAndRecipTerms();
  SystemPotential pot = system->calcEnergy.SystemTotal();

  std::cout
      << "================================================================="
      << std::endl << "-------------------------" << std::endl
      << " STEP: " << step + 1
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
      << "================================================================"
      << std::endl << std::endl;

}

#endif
