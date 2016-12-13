#include "Forcefield.h" //Header spec.
//Setup partner classes
#include "Setup.h"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>

const double BOLTZMANN = 0.0019872041;

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
  if( angles!= NULL)
    delete angles;

}

void Forcefield::Init(const Setup& set)
{
  InitBasicVals(set.config.sys, set.config.in.ffKind);
  particles->Init(set.ff.mie, set.ff.nbfix, set.config.sys,
                  set.config.in.ffKind);
  bonds.Init(set.ff.bond);
  angles->Init(set.ff.angle);
  dihedrals.Init(set.ff.dih);
}

void Forcefield::InitBasicVals(config_setup::SystemVals const& val,
                               config_setup::FFKind const& ffKind)
{
  useLRC = val.ff.doTailCorr;
  T_in_K = val.T.inKelvin;
  rCut = val.ff.cutoff;
  rCutSq = rCut * rCut;
  rCutOver2 = rCut / 2.0;
  scl_14 = val.ff.oneFourScale;
  beta = 1/T_in_K;

  vdwKind = val.ff.VDW_KIND;
  exckind = val.exclude.EXCLUDE_KIND;

  electrostatic = val.elect.enable;
  ewald = val.elect.ewald;
  alpha= val.elect.alpha;
  recip_rcut = val.elect.recip_rcut;


  if(vdwKind == val.ff.VDW_STD_KIND)
    particles = new FFParticle();
  else if(vdwKind == val.ff.VDW_SHIFT_KIND)
    particles = new FF_SHIFT();
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && ffKind.isMARTINI)
    particles = new FF_SWITCH_MARTINI();
  else
    particles = new FF_SWITCH();

  if(ffKind.isMARTINI)
    angles = new FFAngleMartini();
  else
    angles = new FFAngles();

  // Define type of interaction to be included. ex. 1-3, 1-4 and more
  if(exckind == val.exclude.EXC_ONETWO_KIND && ffKind.isMARTINI)
  {
    OneThree = true, OneFour = true, OneN = true;
    std::cout <<
              "REMINDER! 1-3 Interaction and more is ON for Martini forcefield\n"
              << std::endl;
  }
  else if(exckind == val.exclude.EXC_ONEFOUR_KIND)
  {
    OneThree = false, OneFour = false, OneN = true;
    std::cout << "REMINDER! 1-3 and 1-4 Interaction is OFF\n" << std::endl;
  }
  else
  {
    OneThree = false, OneFour = true, OneN = true;
    std::cout << "REMINDER! 1-4 Interaction and more is ON\n" << std::endl;
  }

}
