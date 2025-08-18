/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "System.h"

#include "CalculateEnergy.h"
#include "ConfigSetup.h" //For types directly read from config. file
#include "CrankShaft.h"
#include "EnergyTypes.h"
#include "EnsemblePreprocessor.h"
#include "Ewald.h"
#include "EwaldCached.h"
#include "GOMCEventsProfile.h"
#include "IntraMoleculeExchange1.h"
#include "IntraMoleculeExchange2.h"
#include "IntraMoleculeExchange3.h"
#include "IntraSwap.h"
#include "IntraTargetedSwap.h"
#include "MoleculeExchange1.h"
#include "MoleculeExchange2.h"
#include "MoleculeExchange2Liq.h"
#include "MoleculeExchange3.h"
#include "MoleculeExchange3Liq.h"
#include "MoleculeTransfer.h"
#include "Molecules.h" //For indexing molecules.
#include "MoveBase.h"  //For move bases....
#include "MoveConst.h" //For array of move objects.
#include "MultiParticle.h"
#include "MultiParticleBrownianMotion.h"
#include "NeMTMC.h"
#include "NoEwald.h"
#include "Regrowth.h"
#include "Rotation.h"
#include "Setup.h" //For source of setup data.
#include "StaticVals.h"
#include "TargetedSwap.h"
#include "Translate.h"
#include "VolumeTransfer.h"

System::System(StaticVals &statics, Setup &set, ulong &startStep,
               MultiSim const *const &multisim)
    : statV(statics), boxDimRef(*BoxDim(statics.isOrthogonal)),
#ifdef VARIABLE_PARTICLE_NUMBER
      molLookupRef(molLookup),
#else
      molLookupRef(statics.molLookup),
#endif
      prng(molLookupRef),
#if GOMC_LIB_MPI
      ms(multisim),
#endif
      moveSettings(boxDimRef), cellList(statics.mol, boxDimRef),
      coordinates(boxDimRef, com, molLookupRef, prng, statics.mol),
      com(boxDimRef, coordinates, molLookupRef, statics.mol),
      calcEnergy(statics, *this),
      checkpointSet(startStep, trueStep, molLookupRef, moveSettings,
                    statics.mol, prng, r123wrapper, set),
      vel(statics.forcefield, molLookupRef, statics.mol, prng),
      restartFromCheckpoint(set.config.in.restart.restartFromCheckpoint),
      startStepRef(startStep), trueStep(0) {
  calcEwald = NULL;
#if GOMC_LIB_MPI
  if (ms->parallelTemperingEnabled)
    prngParallelTemp = new PRNG(molLookupRef);
#endif
  prng.Init(set.prng.prngMaker.prng);
  r123wrapper.SetRandomSeed(set.config.in.prng.seed);
#if GOMC_LIB_MPI
  if (ms->parallelTemperingEnabled)
    prngParallelTemp->Init(set.prngParallelTemp.prngMaker.prng);
#endif
  // At this point see if checkpoint is enabled. if so re-initialize
  // step, movesettings, prng, original molecule start and kindex arrays, and
  // original molecule trajectory indices
  if (restartFromCheckpoint) {
    checkpointSet.loadCheckpointFile();
  }
  // overwrite startStepRef with initStep
  if (set.config.sys.step.initStepRead) {
    startStepRef = set.config.sys.step.initStep;
  }
}

System::~System() {
  if (boxDimensions != NULL)
    delete boxDimensions;
  if (calcEwald != NULL)
    delete calcEwald;
  for (int m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
    delete moves[m];
  }
#if GOMC_LIB_MPI
  if (ms->parallelTemperingEnabled)
    delete prngParallelTemp;
#endif
}

void System::Init(Setup &set) {
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(statV.mol, set.pdb.atoms, statV.forcefield,
                 set.config.in.restart.restartFromCheckpoint);
#endif
  moveSettings.Init(statV, set.pdb.remarks, molLookupRef.GetNumKind(),
                    set.config.in.restart.restartFromCheckpoint);
  // allocate memory for atom's velocity if we read the binVelocities
  vel.Init(set.pdb.atoms, set.config.in);
  GOMC_EVENT_START(1, GomcProfileEvent::READ_INPUT_FILES);
  // set coordinates and velocities for atoms in system
  xsc.Init(set.pdb, vel, set.config.in, molLookupRef, statV.mol);
  GOMC_EVENT_STOP(1, GomcProfileEvent::READ_INPUT_FILES);
  boxDimensions->Init(set.config.in.restart, set.config.sys.volume,
                      set.pdb.cryst, statV.forcefield);
  coordinates.InitFromPDB(set.pdb.atoms);
  com.CalcCOM();
  // Allocate space for atom forces
  atomForceRef.Init(set.pdb.atoms.beta.size());
  molForceRef.Init(com.Count());
  // Allocate space for reciprocal force
  atomForceRecRef.Init(set.pdb.atoms.beta.size());
  molForceRecRef.Init(com.Count());
  cellList.SetCutoff();
  cellList.GridAll(boxDimRef, coordinates, molLookupRef);

  // check if we have to use cached version of Ewald or not.
  bool ewald = set.config.sys.elect.ewald;

#ifdef GOMC_CUDA
  if (ewald)
    calcEwald = new Ewald(statV, *this);
  else
    calcEwald = new NoEwald(statV, *this);
#else
  bool cached = set.config.sys.elect.cache;
  if (ewald && cached)
    calcEwald = new EwaldCached(statV, *this);
  else if (ewald && !cached)
    calcEwald = new Ewald(statV, *this);
  else
    calcEwald = new NoEwald(statV, *this);
#endif

  // Initialize lambda before calling SystemTotal
  InitLambda();
  calcEnergy.Init(*this);
  calcEwald->Init();
  potential = calcEnergy.SystemTotal();
  InitMoves(set);
  for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++)
    moveTime[m] = 0.0;
}

void System::InitOver(Setup &set, Molecules &molRef) {
  if (restartFromCheckpoint)
    checkpointSet.InitOver();
}

void System::InitMoves(Setup const &set) {
  moves[mv::DISPLACE] = new Translate(*this, statV);
  moves[mv::MULTIPARTICLE] = new MultiParticle(*this, statV);
  moves[mv::MULTIPARTICLE_BM] = new MultiParticleBrownian(*this, statV);
  moves[mv::ROTATE] = new Rotate(*this, statV);
  moves[mv::INTRA_SWAP] = new IntraSwap(*this, statV);
  moves[mv::REGROWTH] = new Regrowth(*this, statV);
  moves[mv::CRANKSHAFT] = new CrankShaft(*this, statV);
  moves[mv::INTRA_TARGETED_SWAP] = new IntraTargetedSwap(*this, statV);
  if (set.config.sys.intraMemcVal.MEMC1) {
    moves[mv::INTRA_MEMC] = new IntraMoleculeExchange1(*this, statV);
  } else if (set.config.sys.intraMemcVal.MEMC2) {
    moves[mv::INTRA_MEMC] = new IntraMoleculeExchange2(*this, statV);
  } else {
    moves[mv::INTRA_MEMC] = new IntraMoleculeExchange3(*this, statV);
  }

#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  moves[mv::VOL_TRANSFER] = new VolumeTransfer(*this, statV);
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  moves[mv::MOL_TRANSFER] = new MoleculeTransfer(*this, statV);
  if (set.config.sys.memcVal.MEMC1) {
    moves[mv::MEMC] = new MoleculeExchange1(*this, statV);
  } else if (set.config.sys.memcVal.MEMC2) {
    moves[mv::MEMC] = new MoleculeExchange2(*this, statV);
  } else if (set.config.sys.memcVal.MEMC2Liq) {
    moves[mv::MEMC] = new MoleculeExchange2Liq(*this, statV);
  } else if (set.config.sys.memcVal.MEMC3Liq) {
    moves[mv::MEMC] = new MoleculeExchange3Liq(*this, statV);
  } else {
    moves[mv::MEMC] = new MoleculeExchange3(*this, statV);
  }
  moves[mv::NE_MTMC] = new NEMTMC(*this, statV);
  moves[mv::TARGETED_SWAP] = new TargetedSwap(*this, statV);

#endif
}

void System::InitLambda() {
  // Set the pointer to cudaVar and initialize the values
  lambdaRef.Init(
#ifdef GOMC_CUDA
      statV.forcefield.particles->getCUDAVars()
#endif
  );

  if (statV.freeEnVal.enable) {
    bool found = false;
    for (uint k = 0; k < statV.mol.GetKindsCount(); k++) {
      std::string kindName = statV.mol.kinds[k].name;
      if (statV.freeEnVal.molType == kindName) {
        found = true;
        uint totalMol = molLookupRef.NumKindInBox(k, mv::BOX0);
        // In PDB file, molIndex start from 1.
        uint FEmolIndex = statV.freeEnVal.molIndex - 1;
        if (totalMol == 0) {
          found = false;
        } else if (totalMol <= FEmolIndex) {
          std::cout << "Error: Molecule index " << statV.freeEnVal.molIndex
                    << " of kind " << kindName
                    << " does not exist in the simulation box!\n";
          exit(EXIT_FAILURE);
        } else {
          uint m = molLookupRef.GetMolNum(FEmolIndex, k, mv::BOX0);
          uint state = statV.freeEnVal.iState;
          double lambdaCoulomb = statV.freeEnVal.lambdaCoulomb[state];
          double lambdaVDW = statV.freeEnVal.lambdaVDW[state];
          lambdaRef.Set(lambdaVDW, lambdaCoulomb, m, k, mv::BOX0);
        }
        break;
      }
    }

    if (!found) {
      std::cout << "Error: No molecule of kind " << statV.freeEnVal.molType
                << " in the simulation box! \n";
      exit(EXIT_FAILURE);
    }
  }
}

void System::RecalculateTrajectory(Setup &set, uint frameNum) {
  set.pdb.Init(set.config.in.restart, set.config.in.files.pdb.name, frameNum);
  statV.InitOver(set, *this);
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(statV.mol, set.pdb.atoms, statV.forcefield,
                 set.config.in.restart.restartFromCheckpoint);
#endif
  coordinates.InitFromPDB(set.pdb.atoms);
  com.CalcCOM();
  cellList.GridAll(boxDimRef, coordinates, molLookupRef);
  calcEnergy.Init(*this);
  calcEwald->Init();
  potential = calcEnergy.SystemTotal();
}

void System::ChooseAndRunMove(const ulong step) {
  if (restartFromCheckpoint) {
    // step - startStepRef = Δstep
    // true step + Δstep
    r123wrapper.SetStep(trueStep + step - startStepRef);
  } else {
    r123wrapper.SetStep(step);
  }
  double draw = 0;
  uint majKind = 0;
  PickMove(majKind, draw);
#ifndef NDEBUG
  std::cout << "Step " << step + 1 << ": Picked move #" << majKind << ": "
            << str::MoveTypetoStr(majKind) << " move" << std::endl;
#endif
  time.SetStart();
  RunMove(majKind, draw, step);
  time.SetStop();
  moveTime[majKind] += time.GetTimDiff();
}
void System::PickMove(uint &kind, double &draw) {
  prng.PickArbDist(kind, draw, statV.movePerc, statV.totalPerc,
                   mv::MOVE_KINDS_TOTAL);
}

void System::RunMove(uint majKind, double draw, const ulong step) {
  // return now if move targets molecule and there's none in that box.
  uint rejectState = SetParams(majKind, draw);
  // If single atom, redo move as displacement
  if (rejectState == mv::fail_state::ROTATE_ON_SINGLE_ATOM) {
    majKind = mv::DISPLACE;
    Translate *disp = static_cast<Translate *>(moves[mv::DISPLACE]);
    Rotate *rot = static_cast<Rotate *>(moves[mv::ROTATE]);
    rejectState = disp->ReplaceRot(*rot);
  }
  if (rejectState == mv::fail_state::NO_FAIL)
    rejectState = Transform(majKind);
  if (rejectState == mv::fail_state::NO_FAIL)
    CalcEn(majKind);
  Accept(majKind, rejectState, step);
}

uint System::SetParams(const uint kind, const double draw) {
  return moves[kind]->Prep(draw, statV.movePerc[kind]);
}

uint System::Transform(const uint kind) { return moves[kind]->Transform(); }

void System::CalcEn(const uint kind) { moves[kind]->CalcEn(); }

void System::Accept(const uint kind, const uint rejectState, const ulong step) {
  moves[kind]->Accept(rejectState, step);
}

void System::PrintAcceptance() {
  std::cout << std::endl;
  printf("%-24s %-15s", "Move Type", "Mol. Kind");
  for (uint b = 0; b < BOX_TOTAL; b++) {
    sstrm::Converter toStr;
    std::string numStr = "";
    toStr << b;
    toStr >> numStr;
    numStr = "BOX_" + numStr;
    printf("%-11s", numStr.c_str());
  }
  std::cout << std::endl;

  for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
    if (statV.movePerc[m] > 0.0)
      moves[m]->PrintAcceptKind();
  }
  std::cout << std::endl;
}

void System::PrintTime() {
  std::string moveName;
  // std::cout << "MC moves Execution time:\n";
  for (auto m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
    moveName = str::MoveTypetoStr(m) + ':';
    printf("%-36s %10.4f    sec.\n", moveName.c_str(), moveTime[m]);
  }
  std::cout << std::endl;
}
