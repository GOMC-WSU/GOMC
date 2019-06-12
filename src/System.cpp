/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "EnsemblePreprocessor.h"
#include "System.h"
#include "CalculateEnergy.h"
#include "EwaldCached.h"
#include "Ewald.h"
#include "NoEwald.h"
#include "EnergyTypes.h"
#include "Setup.h"               //For source of setup data.
#include "ConfigSetup.h"         //For types directly read from config. file
#include "StaticVals.h"
#include "Molecules.h"           //For indexing molecules.
#include "MoveConst.h"           //For array of move objects.
#include "MoveBase.h"            //For move bases....
#include "Rotation.h"
#include "Translate.h"
#include "VolumeTransfer.h"
#include "MoleculeTransfer.h"
#include "IntraSwap.h"
#include "MultiParticle.h"
#include "Regrowth.h"
#include "MoleculeExchange1.h"
#include "MoleculeExchange2.h"
#include "MoleculeExchange3.h"
#include "IntraMoleculeExchange1.h"
#include "IntraMoleculeExchange2.h"
#include "IntraMoleculeExchange3.h"
#include "CrankShaft.h"

System::System(StaticVals& statics) :
  statV(statics),
#ifdef VARIABLE_VOLUME
  boxDimRef(*BoxDim(statics.isOrthogonal)),
#else
  boxDimRef(*statics.GetBoxDim()),
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookupRef(molLookup),
#else
  molLookupRef(statics.molLookup),
#endif
  prng(molLookupRef),
  coordinates(boxDimRef, com, molLookupRef, prng, statics.mol),
  com(boxDimRef, coordinates, molLookupRef, statics.mol),
  moveSettings(boxDimRef), cellList(statics.mol, boxDimRef),
  calcEnergy(statics, *this), checkpointSet(*this, statics)
{
  calcEwald = NULL;
}

System::~System()
{
#ifdef VARIABLE_VOLUME
  if (boxDimensions != NULL)
    delete boxDimensions;
#endif
  if (calcEwald != NULL)
    delete calcEwald;
  delete moves[mv::DISPLACE];
  delete moves[mv::ROTATE];
  delete moves[mv::INTRA_SWAP];
  delete moves[mv::REGROWTH];
  delete moves[mv::INTRA_MEMC];
  delete moves[mv::CRANKSHAFT];
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  delete moves[mv::VOL_TRANSFER];
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  delete moves[mv::MOL_TRANSFER];
  delete moves[mv::MEMC];
#endif
}

void System::Init(Setup const& set, ulong & startStep)
{
  prng.Init(set.prng.prngMaker.prng);
#ifdef VARIABLE_VOLUME
  boxDimensions->Init(set.config.in.restart,
                      set.config.sys.volume, set.pdb.cryst,
                      statV.forcefield);
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(statV.mol, set.pdb.atoms);
#endif
  moveSettings.Init(statV, set.pdb.remarks, molLookupRef.GetNumKind());
  //Note... the following calls use box iterators, so must come after
  //the molecule lookup initialization, in case we're in a constant
  //particle/molecule ensemble, e.g. NVT
  coordinates.InitFromPDB(set.pdb.atoms);

  // At this point see if checkpoint is enabled. if so re-initialize
  // coordinates, prng, mollookup, step, boxdim, and movesettings
  if(set.config.in.restart.restartFromCheckpoint) {
    checkpointSet.ReadAll();
    checkpointSet.SetStepNumber(startStep);
    checkpointSet.SetBoxDimensions(boxDimRef);
    checkpointSet.SetPRNGVariables(prng);
    checkpointSet.SetCoordinates(coordinates);
    checkpointSet.SetMoleculeLookup(molLookupRef);
    checkpointSet.SetMoveSettings(moveSettings);
  }

  com.CalcCOM();
  // Allocate space for atom forces
  atomForceRef.Init(set.pdb.atoms.beta.size());
  molForceRef.Init(com.Count());
  // Allocate space for reciprocate force
  atomForceRecRef.Init(set.pdb.atoms.beta.size());
  molForceRecRef.Init(com.Count());
  cellList.SetCutoff();
  cellList.GridAll(boxDimRef, coordinates, molLookupRef);

  //check if we have to use cached version of ewlad or not.
  bool ewald = set.config.sys.elect.ewald;
  bool cached = set.config.sys.elect.cache;

#ifdef GOMC_CUDA
  if(ewald)
    calcEwald = new Ewald(statV, *this);
  else
    calcEwald = new NoEwald(statV, *this);
#else
  if (ewald && cached)
    calcEwald = new EwaldCached(statV, *this);
  else if (ewald && !cached)
    calcEwald = new Ewald(statV, *this);
  else
    calcEwald = new NoEwald(statV, *this);
#endif

  calcEnergy.Init(*this);
  calcEwald->Init();
  potential = calcEnergy.SystemTotal();
  InitMoves(set);
  for(uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++)
    moveTime[m] = 0.0;
}

void System::InitMoves(Setup const& set)
{
  moves[mv::DISPLACE] = new Translate(*this, statV);
  moves[mv::MULTIPARTICLE] = new MultiParticle(*this, statV);
  moves[mv::ROTATE] = new Rotate(*this, statV);
  moves[mv::INTRA_SWAP] = new IntraSwap(*this, statV);
  moves[mv::REGROWTH] = new Regrowth(*this, statV);
  moves[mv::CRANKSHAFT] = new CrankShaft(*this, statV);
  if(set.config.sys.intraMemcVal.MEMC1) {
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
  if(set.config.sys.memcVal.MEMC1) {
    moves[mv::MEMC] = new MoleculeExchange1(*this, statV);
  } else if (set.config.sys.memcVal.MEMC2) {
    moves[mv::MEMC] = new MoleculeExchange2(*this, statV);
  } else {
    moves[mv::MEMC] = new MoleculeExchange3(*this, statV);
  }
#endif
}

void System::RecalculateTrajectory(Setup &set, uint frameNum)
{
  set.pdb.Init(set.config.in.restart, set.config.in.files.pdb.name, frameNum);
  statV.InitOver(set, *this);
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(statV.mol, set.pdb.atoms);
#endif
  coordinates.InitFromPDB(set.pdb.atoms);
  com.CalcCOM();
  cellList.GridAll(boxDimRef, coordinates, molLookupRef);
  calcEnergy.Init(*this);
  calcEwald->Init();
  potential = calcEnergy.SystemTotal();
}

void System::ChooseAndRunMove(const uint step)
{
  real draw = 0;
  uint majKind = 0;
  PickMove(majKind, draw);
  time.SetStart();
  RunMove(majKind, draw, step);
  time.SetStop();
  moveTime[majKind] += time.GetTimDiff();
}
void System::PickMove(uint & kind, real & draw)
{
  prng.PickArbDist(kind, draw, statV.movePerc, statV.totalPerc,
                   mv::MOVE_KINDS_TOTAL);
}

void System::RunMove(uint majKind, real draw, const uint step)
{
  //return now if move targets molecule and there's none in that box.
  uint rejectState = SetParams(majKind, draw);
  //If single atom, redo move as displacement
  if (rejectState == mv::fail_state::ROTATE_ON_SINGLE_ATOM) {
    majKind = mv::DISPLACE;
    Translate * disp = static_cast<Translate *>(moves[mv::DISPLACE]);
    Rotate * rot = static_cast<Rotate *>(moves[mv::ROTATE]);
    rejectState = disp->ReplaceRot(*rot);
  }
  if (rejectState == mv::fail_state::NO_FAIL)
    rejectState = Transform(majKind);
  if (rejectState == mv::fail_state::NO_FAIL)
    CalcEn(majKind);
  Accept(majKind, rejectState, step);
}

uint System::SetParams(const uint kind, const real draw)
{
  return moves[kind]->Prep(draw, statV.movePerc[kind]);
}

uint System::Transform(const uint kind)
{
  return moves[kind]->Transform();
}

void System::CalcEn(const uint kind)
{
  moves[kind]->CalcEn();
}

void System::Accept(const uint kind, const uint rejectState, const uint step)
{
  moves[kind]->Accept(rejectState, step);
}

void System::PrintAcceptance()
{
  std::cout << std::endl;
  printf("%-24s %-15s", "Move Type", "Mol. Kind");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    sstrm::Converter toStr;
    std::string numStr = "";
    toStr << b;
    toStr >> numStr;
    numStr = "BOX_" + numStr;
    printf("%-11s", numStr.c_str());
  }
  std::cout << std::endl;

  for(uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
    if(statV.movePerc[m] > 0.0)
      moves[m]->PrintAcceptKind();
  }
  std::cout << std::endl;
}

void System::PrintTime()
{
  //std::cout << "MC moves Execution time:\n";
  printf("%-36s %10.4f    sec.\n", "Displacement:", moveTime[mv::DISPLACE]);
  printf("%-36s %10.4f    sec.\n", "Rotation:", moveTime[mv::ROTATE]);
  printf("%-36s %10.4f    sec.\n", "MultiParticle:", moveTime[mv::MULTIPARTICLE]);
  printf("%-36s %10.4f    sec.\n", "Intra-Swap:", moveTime[mv::INTRA_SWAP]);
  printf("%-36s %10.4f    sec.\n", "Regrowth:", moveTime[mv::REGROWTH]);
  printf("%-36s %10.4f    sec.\n", "Intra-MEMC:", moveTime[mv::INTRA_MEMC]);
  printf("%-36s %10.4f    sec.\n", "Crank-Shaft:", moveTime[mv::CRANKSHAFT]);

#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  printf("%-36s %10.4f    sec.\n", "Mol-Transfer:",
         moveTime[mv::MOL_TRANSFER]);
  printf("%-36s %10.4f    sec.\n", "MEMC:", moveTime[mv::MEMC]);
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  printf("%-36s %10.4f    sec.\n", "Vol-Transfer:", moveTime[mv::VOL_TRANSFER]);
#endif
  std::cout << std::endl;
}
