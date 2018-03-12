/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
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
#include "MoleculeTransfer.h"
#include "IntraSwap.h"

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
  calcEnergy(statics, *this)
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
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  delete moves[mv::VOL_TRANSFER];
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  delete moves[mv::MOL_TRANSFER];
#endif
}

void System::Init(Setup const& set)
{
  prng.Init(set.prng.prngMaker.prng);
#ifdef VARIABLE_VOLUME
  boxDimensions->Init(set.config.in.restart,
                      set.config.sys.volume, set.pdb.cryst,
                      statV.forcefield.rCut,
                      statV.forcefield.rCutSq);
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
  molLookup.Init(statV.mol, set.pdb.atoms);
#endif
  moveSettings.Init(statV, set.pdb.remarks);
  //Note... the following calls use box iterators, so must come after
  //the molecule lookup initialization, in case we're in a constant
  //particle/molecule ensemble, e.g. NVT
  coordinates.InitFromPDB(set.pdb.atoms);
  com.CalcCOM();
  cellList.SetCutoff(statV.forcefield.rCut);
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
  InitMoves();
}

void System::InitMoves()
{
  moves[mv::DISPLACE] = new Translate(*this, statV);
  moves[mv::ROTATE] = new Rotate(*this, statV);
  moves[mv::INTRA_SWAP] = new IntraSwap(*this, statV);
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
  moves[mv::VOL_TRANSFER] = new VolumeTransfer(*this, statV);
#endif
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
  moves[mv::MOL_TRANSFER] = new MoleculeTransfer(*this, statV);
#endif
}

void System::ChooseAndRunMove(const uint step)
{
  double draw = 0;
  uint majKind = 0;
  PickMove(majKind, draw);
  RunMove(majKind, draw, step);
}
void System::PickMove(uint & kind, double & draw)
{
  prng.PickArbDist(kind, draw, statV.movePerc, statV.totalPerc,
                   mv::MOVE_KINDS_TOTAL);
}

void System::RunMove(uint majKind, double draw, const uint step)
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

uint System::SetParams(const uint kind, const double draw)
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
