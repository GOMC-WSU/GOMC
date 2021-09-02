/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Simulation.h"
#include "Setup.h"          //For setup object
#include "NumLib.h"
#include "EnergyTypes.h"
#include <iostream>
#include <iomanip>
#include "CUDAMemoryManager.cuh"
#include "GOMCEventsProfile.h"

#define EPSILON 0.001

Simulation::Simulation(char const*const configFileName, MultiSim const*const& multisim): ms(multisim)
{
  GOMC_EVENT_START(1, GomcProfileEvent::INITIALIZE);
  startStep = 0;
  GOMC_EVENT_START(1, GomcProfileEvent::READ_INPUT_FILES);
  set.Init(configFileName, multisim);
  GOMC_EVENT_STOP(1, GomcProfileEvent::READ_INPUT_FILES);
  totalSteps = set.config.sys.step.total;
  staticValues = new StaticVals(set);
  system = new System(*staticValues, set, multisim);
  staticValues->Init(set, *system);
  system->Init(set, startStep);
  //recalc Init for static value for initializing ewald since ewald is
  //initialized in system
  staticValues->InitOver(set, *system);
  system->InitOver(set, staticValues->mol);
  cpu = new CPUSide(*system, *staticValues, set);
  cpu->Init(set.pdb, set.config.out, set.config.sys.step.equil,
            totalSteps, startStep);
            
  if(totalSteps == 0) {
    frameSteps = set.pdb.GetFrameSteps(set.config.in.files.pdb.name);
  }
#if GOMC_LIB_MPI
  // set.config.sys.step.parallelTemp is a boolean for enabling/disabling parallel tempering
  PTUtils = set.config.sys.step.parallelTemp ? new ParallelTemperingUtilities(ms, *system, *staticValues, set.config.sys.step.parallelTempFreq, set.config.sys.step.parallelTemperingAttemptsPerExchange) : NULL;
  exchangeResults.resize(ms->worldSize, false);
#endif
  GOMC_EVENT_STOP(1, GomcProfileEvent::INITIALIZE);
}

Simulation::~Simulation()
{
  GOMC_EVENT_START(1, GomcProfileEvent::DESTRUCTION);
  delete cpu;
  delete system;
  delete staticValues;
#ifdef GOMC_CUDA
  CUDAMemoryManager::isFreed();
#endif
  GOMC_EVENT_STOP(1, GomcProfileEvent::DESTRUCTION);
}

void Simulation::RunSimulation(void)
{
  GOMC_EVENT_START(1, GomcProfileEvent::MC_RUN);
  double startEnergy = system->potential.totalEnergy.total;
  if(totalSteps == 0) {
    for(int i = 0; i < (int) frameSteps.size(); i++) {
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
      if(std::abs(currEnergy - startEnergy) > 1.0e+10) {
        printf("Info: Recalculating the total energies to insure the accuracy"
               " of the computed \n"
               "      running energies.\n\n");
        system->calcEwald->UpdateVectorsAndRecipTerms(true);
        system->potential = system->calcEnergy.SystemTotal();
      }
    }

#if GOMC_LIB_MPI
    //
    if(staticValues->simEventFreq.parallelTemp && step > cpu->equilSteps && step % staticValues->simEventFreq.parallelTempFreq == 0) {

      int maxSwap = 0;
      /* Number of rounds of exchanges needed to deal with any multiple
      * exchanges. */
      /* Where each replica ends up after the exchange attempt(s). */
      /* The order in which multiple exchanges will occur. */
      bool bThisReplicaExchanged = false;

      system->potential = system->calcEnergy.SystemTotal();
      PTUtils->evaluateExchangeCriteria(step);
      PTUtils->prepareToDoExchange(ms->worldRank, &maxSwap, &bThisReplicaExchanged);
      PTUtils->conductExchanges(ms->worldRank, system->coordinates, system->com, maxSwap, bThisReplicaExchanged);
    }

#endif

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RecalculateAndCheck();
#endif
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::MC_RUN);
  if(!RecalculateAndCheck()) {
    std::cerr << "Warning: Updated energy differs from Recalculated Energy!\n";
  }
  system->PrintAcceptance();
  system->PrintTime();
#if GOMC_LIB_MPI
  if (staticValues->simEventFreq.parallelTemp)
    PTUtils->print_replica_exchange_statistics(ms->fplog);
#endif
}

bool Simulation::RecalculateAndCheck(void)
{
  system->calcEwald->UpdateVectorsAndRecipTerms(false);
  SystemPotential pot = system->calcEnergy.SystemTotal();

  bool compare = true;
  compare &= num::approximatelyEqual(system->potential.totalEnergy.intraBond, pot.totalEnergy.intraBond, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.intraNonbond, pot.totalEnergy.intraNonbond, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.inter, pot.totalEnergy.inter, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.tc, pot.totalEnergy.tc, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.real, pot.totalEnergy.real, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.self, pot.totalEnergy.self, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.correction, pot.totalEnergy.correction, EPSILON);
  compare &= num::approximatelyEqual(system->potential.totalEnergy.recip, pot.totalEnergy.recip, EPSILON);

  if(!compare) {
    std::cout
        << "=================================================================\n"
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

  return compare;
}
  #if GOMC_GTEST_MPI || GOMC_GTEST
      double Simulation::GetSystemEnergy(void){
        return system->potential.Total();
      }
      Coordinates& Simulation::GetCoordinates(void){
        return system->coordinates;
      }
      COM& Simulation::GetCOMs(void){
        return system->com;
      }
      CellList& Simulation::GetCellList(void){
        return system->cellList;
      }
      StaticVals* Simulation::GetStaticVals(void){
        return staticValues;
      }
      CalculateEnergy& Simulation::GetCalcEn(){
        return system->calcEnergy;
      }

  #endif

  #if GOMC_GTEST_MPI
    #if ENSEMBLE == NPT
        void Simulation::SetGlobalVolumes(int worldRank){
          std::cout << "I (" << worldRank << ") think my volume is " 
            << GetVolume(0) << std::endl;
          PTUtils->SetGlobalVolumes(worldRank, GetVolume(0));
           std::cout << "My (" << worldRank << ") volume is now " 
            << GetVolume(0) << std::endl;
        }

        double Simulation::GetVolume(int box){
          return system->boxDimRef.volume[box];
        }
    #elif ENSEMBLE == GCMC
        void Simulation::SetGlobalChemPots(int worldRank){
          std::cout << "I (" << worldRank << ") think my chem pots are " 
            << GetChemPot(0) << std::endl;
          PTUtils->SetGlobalChemPots(worldRank, *system, *staticValues);
           std::cout << "My (" << worldRank << ") volume is now " 
            << GetChemPot(0) << std::endl;
        }

        double Simulation::GetChemPot(int kind){
          return staticValues->mol.kinds[kind].chemPot;
        }
        void Simulation::SetGlobalNumberOfMolecules(int worldRank){
          std::cout << "I (" << worldRank << ") think my molecule numbers are " 
            << GetChemPot(0) << std::endl;
          PTUtils->SetGlobalNumberOfMolecules(worldRank, *system, *staticValues);
           std::cout << "My (" << worldRank << ") volume is now " 
            << GetChemPot(0) << std::endl;
        }

        void Simulation::GetNumberOfMolecules(std::vector<uint> & moleculeNumbers){
          for (uint box = 0; box < BOX_TOTAL; ++box)
            for (int kind = 0; kind < staticValues->mol.kindsCount; ++kind)
              moleculeNumbers.push_back(system->molLookupRef.NumKindInBox(kind, box));    
        }
    #endif
      void Simulation::ExchangeReplicas(int worldRank){
        std::cout << "calling fe" << std::endl;
        PTUtils->forceExchange(worldRank, *system, *staticValues);
        std::cout << "returned from fe" << std::endl;
      }
  #endif