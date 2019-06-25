/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Forcefield.h" //Header spec.
//Setup partner classes
#include "Setup.h"
#include "FFShift.h"
#include "FFSwitch.h"
#include "FFSwitchMartini.h"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

Forcefield::Forcefield()
{
  particles = NULL;
  angles = NULL;
  OneThree = false; //default behavior is to turn off 1-3 interaction
  OneFour = true;   // to turn on 1-4 interaction
  OneN = true;      // and turn on 1-n interaction
}

Forcefield::~Forcefield()
{
  if(particles != NULL)
    delete particles;
  if( angles != NULL)
    delete angles;

}

void Forcefield::Init(const Setup& set)
{
  InitBasicVals(set.config.sys, set.config.in.ffKind);
  particles->Init(set.ff.mie, set.ff.nbfix);
  bonds.Init(set.ff.bond);
  angles->Init(set.ff.angle);
  dihedrals.Init(set.ff.dih);
}

void Forcefield::InitBasicVals(config_setup::SystemVals const& val,
                               config_setup::FFKind const& ffKind)
{
  useLRC = val.ff.doTailCorr;
  energyTable = val.ff.energyTable;
  T_in_K = val.T.inKelvin;
  rCut = val.ff.cutoff;
  rCutSq = rCut * rCut;
  rCutLow = val.ff.cutoffLow;
  rCutLowSq = rCutLow * rCutLow;
  scaling_14 = val.elect.oneFourScale;
  beta = 1 / T_in_K;

  vdwKind = val.ff.VDW_KIND;
  exckind = val.exclude.EXCLUDE_KIND;

  electrostatic = val.elect.enable;
  ewald = val.elect.ewald;
  tolerance = val.elect.tolerance;
  rswitch = val.ff.rswitch;
  dielectric = val.elect.dielectric;

  for(uint b = 0 ; b < BOX_TOTAL; b++) {
    rCutCoulomb[b] = val.elect.cutoffCoulomb[b];
    rCutCoulombSq[b] = rCutCoulomb[b] * rCutCoulomb[b];
    alpha[b] = sqrt(-log(tolerance)) / rCutCoulomb[b];
    alphaSq[b] = alpha[b] * alpha[b];
    recip_rcut[b] = -2.0 * log(tolerance) / rCutCoulomb[b];
    recip_rcut_Sq[b] = recip_rcut[b] * recip_rcut[b];
  }

  vdwGeometricSigma = val.ff.vdwGeometricSigma;
  isMartini = ffKind.isMARTINI;

#if ENSEMBLE == GCMC
  isFugacity = val.chemPot.isFugacity;
#endif

  if(vdwKind == val.ff.VDW_STD_KIND)
    particles = new FFParticle(*this);
  else if(vdwKind == val.ff.VDW_SHIFT_KIND)
    particles = new FF_SHIFT(*this);
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && ffKind.isMARTINI)
    particles = new FF_SWITCH_MARTINI(*this);
  else
    particles = new FF_SWITCH(*this);

  if(ffKind.isMARTINI)
    angles = new FFAngleMartini();
  else
    angles = new FFAngles();

  // Define type of interaction to be included. ex. 1-3, 1-4 and more
  if(exckind == val.exclude.EXC_ONETWO_KIND) {
    OneThree = true, OneFour = true, OneN = true;
  } else if(exckind == val.exclude.EXC_ONETHREE_KIND) {
    OneThree = false, OneFour = true, OneN = true;
  } else if(exckind == val.exclude.EXC_ONEFOUR_KIND) {
    OneThree = false, OneFour = false, OneN = true;
  } else {
    std::cout << "Error: Unknown exclude value.\n";
    exit(EXIT_FAILURE);
  }

}
