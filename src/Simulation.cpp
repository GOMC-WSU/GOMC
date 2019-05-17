/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
#include "ReplDirSetup.h"

Simulation::Simulation(char const*const configFileName)
{
  startStep = 0;
  this->configFileName = configFileName;
  //NOTE:
  //IMPORTANT! Keep this order...
  //as system depends on staticValues, and cpu sometimes depends on both.
  set.Init(configFileName);
  totalSteps = set.config.sys.step.total;
  absoluteTotalSteps = set.config.sys.step.total;
  staticValues = new StaticVals(set);
  system = new System(*staticValues);
  staticValues->Init(set, *system);
  system->Init(set, startStep);
  initReplExParams(set.config.in.replValues);
  //recal Init for static value for initializing ewald since ewald is
  //initialized in system
  staticValues->InitOver(set, *system);

  /* "Replica_Exchange_Simulation" is the default value, so if compare returns a value !=0 it means the user
      provided a title for the multisimulation.  Or, the user provided an exchange interval > 0 which implies
      a multisimulation.  The title rule allows for multisimulations to be run without exchanging, which will be
      useful for testing the same system with/without exchanging.
  */
  if (replExParams.multiSimTitle.compare("Replica_Exchange_Simulation") || replExParams.exchangeInterval > 0){
    setupHierarchicalDirectoryStructure();
  }

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
}

void Simulation::RunSimulation(void)
{
  if (startStep == 0)
    startEnergy = system->potential.totalEnergy.total;
  
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

#ifndef NDEBUG
    if((step + 1) % 1000 == 0)
      RunningCheck(step);
#endif
  }
  if(totalSteps == absoluteTotalSteps){
    if (cpu->getConsoleOutput()->getConsoleToFile() != NULL){
      system->PrintAcceptance(cpu->getConsoleOutput()->getConsoleToFile());
      system->PrintTime(cpu->getConsoleOutput()->getConsoleToFile());
    } else {
      system->PrintAcceptance();
      system->PrintTime();   
    }
  }
}

void Simulation::RunNSteps(ulong NSteps)
{
  ulong temporaryStorageOfTotalSteps = totalSteps;
  totalSteps = NSteps + startStep;
  RunSimulation();
  startStep += NSteps;
  totalSteps = temporaryStorageOfTotalSteps;
}

ulong Simulation::getTotalSteps(){
  return totalSteps;
}

ulong Simulation::getStartStep(){
  return startStep;
}

ulong Simulation::getEquilSteps(){
  return cpu->equilSteps;
}

double Simulation::getT_in_K(){
  return staticValues->forcefield.T_in_K;
}

double Simulation::getBeta(){
  return staticValues->forcefield.beta;
}

ulong Simulation::getExchangeInterval(){
  return replExParams.exchangeInterval;
}

double Simulation::getEpot(){
  return system->potential.totalEnergy.total;
}

double Simulation::getEpotBox(uint i){
  return system->potential.boxEnergy[i].total;
}

CPUSide* Simulation::getCPUSide(){
  return cpu;
}

System* Simulation::getSystem(){
  return system;
}

StaticVals* Simulation::getStaticValues(){
  return staticValues;
}

Clock* Simulation::getClock(){
  return cpu->getClock();
}

std::string Simulation::getConfigFileName(){
  return configFileName;
}

std::string Simulation::getMultiSimTitle(){
  return replExParams.multiSimTitle;
}

int Simulation::getReplExSeed(){
  return replExParams.randomSeed;
}

#if ENSEMBLE == NPT || ENSEMBLE == GEMC

uint Simulation::getKindOfGEMC(){
  return staticValues->kindOfGEMC;
}

double Simulation::getPressure(){
  return staticValues->pressure;
}

double Simulation::getVolume(){
  return *system->boxDimRef.volume; 
}

double Simulation::getVolume(uint i){
  return system->boxDimRef.volume[i];
}
#endif



#if ENSEMBLE == GCMC
int Simulation::getNumOfParticles(uint i){
  uint numMolecules = 0;
  uint box = 0;
  numMolecules = system->molLookup.NumKindInBox(i, box);
  return (int)numMolecules;
}

double Simulation::getChemicalPotential(uint kind){
  double chemPot = 0.0;
  uint box = 0;
  chemPot = staticValues->mol.kinds[kind].chemPot;    
  //chemPot = system->molLookup.NumKindInBox(kind, box) * staticValues->mol.kinds[kind].chemPot;    

  return chemPot;
}
#endif


void Simulation::setT_in_K(double T_in_K){
  staticValues->forcefield.T_in_K = T_in_K;
}

void Simulation::setBeta(double beta){
  staticValues->forcefield.beta = beta;
}

void Simulation::setCPUSide(CPUSide* cpuToSet){
  cpu = cpuToSet;
}

void Simulation::initReplExParams(struct config_setup::ReplicaExchangeValuesFromConf replExValuesFromConfFile){
  replExParams.exchangeInterval = replExValuesFromConfFile.exchangeInterval;
  replExParams.numExchanges = replExValuesFromConfFile.numExchanges;
  replExParams.randomSeed = replExValuesFromConfFile.randomSeed;
  replExParams.multiSimTitle = replExValuesFromConfFile.multiSimTitle;
}

void Simulation::setupHierarchicalDirectoryStructure(){
  ReplDirSetup rd(staticValues->forcefield.T_in_K, replExParams);    
  set.config.out.replica_path =  rd.path_to_replica_directory;    
  set.config.out.useMultidir =  true;
      
  std::stringstream replica_stream;    
  replica_stream << set.config.out.replica_path << set.config.out.state.files.psf.name;    
  set.config.out.state.files.psf.name = replica_stream.str();    

  for(int i = 0; i < BOX_TOTAL; i++) {    
    std::stringstream replica_stream;    
    replica_stream << set.config.out.replica_path << set.config.out.state.files.pdb.name[i];    
    set.config.out.state.files.pdb.name[i] = replica_stream.str();    
  }    
  std::stringstream replica_stream1;    
  replica_stream1 << set.config.out.replica_path << set.config.out.state.files.seed.name;    
  set.config.out.state.files.seed.name = replica_stream1.str();

  std::stringstream replica_stream2;
  replica_stream2 << set.config.out.replica_path << set.config.out.state.files.console.name;
  set.config.out.state.files.console.name = replica_stream2.str();

}

void Simulation::attachNewCPUSideToLocalSysAndStatV(){
  cpu->reInitVarRef(system, staticValues);
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
