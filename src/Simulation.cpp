/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
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
#include "CUDAMemoryManager.cuh"

#define EPSILON 0.001

Simulation::Simulation(char const*const configFileName, MultiSim const*const& multisim)
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
}

Simulation::~Simulation()
{
  delete cpu;
  delete system;
  delete staticValues;
#ifdef GOMC_CUDA
  CUDAMemoryManager::isFreed();
#endif
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
      if(std::abs(currEnergy - startEnergy) > 1.0e+10) {
        printf("Info: Recalculating the total energies to insure the accuracy"
               " of the computed \n"
               "      running energies.\n\n");
        system->calcEwald->Init();
        system->potential = system->calcEnergy.SystemTotal();
      }
    }

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  if(!RecalculateAndCheck()) {
    std::cerr << "Warning: Updated energy differs from Recalculated Energy!\n";
  }
  system->PrintAcceptance();
  system->PrintTime();
}

bool Simulation::RecalculateAndCheck(void)
{
  system->calcEwald->UpdateVectorsAndRecipTerms(false);
  SystemPotential pot = system->calcEnergy.SystemTotal();

  bool compare = true;
  compare &= std::abs(system->potential.totalEnergy.intraBond - pot.totalEnergy.intraBond) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.intraNonbond - pot.totalEnergy.intraNonbond) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.inter - pot.totalEnergy.inter) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.tc - pot.totalEnergy.tc) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.real - pot.totalEnergy.real) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.self - pot.totalEnergy.self) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.correction - pot.totalEnergy.correction) < EPSILON;
  compare &= std::abs(system->potential.totalEnergy.recip - pot.totalEnergy.recip) < EPSILON;

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

#ifndef NDEBUG
void Simulation::RunningCheck(const uint step)
{
  system->calcEwald->UpdateVectorsAndRecipTerms(false);
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
