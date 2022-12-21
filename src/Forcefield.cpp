/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "Forcefield.h" //Header spec.
// Setup partner classes
#include "FFExp6.h"
#include "FFShift.h"
#include "FFSwitch.h"
#include "FFSwitchMartini.h"
#include "Setup.h"
#define _USE_MATH_DEFINES
#include <cmath>

Forcefield::Forcefield() {
  particles = NULL;
  angles = NULL;
  OneThree = false; // default behavior is to turn off 1-3 interaction
  OneFour = true;   // to turn on 1-4 interaction
  OneN = true;      // and turn on 1-n interaction
}

Forcefield::~Forcefield() {
  if (particles != NULL)
    delete particles;
  if (angles != NULL)
    delete angles;
}

void Forcefield::Init(const Setup &set) {
  InitBasicVals(set.config.sys, set.config.in.ffKind);
  particles->Init(set.ff.mie, set.ff.nbfix);
  bonds.Init(set.ff.bond);
  angles->Init(set.ff.angle);
  dihedrals.Init(set.ff.dih);
}

void Forcefield::InitBasicVals(config_setup::SystemVals const &val,
                               config_setup::FFKind const &ffKind) {
  useLRC = val.ff.doTailCorr;
  useIPC = val.ff.doImpulsePressureCorr;
  T_in_K = val.T.inKelvin;
  rCut = val.ff.cutoff;
  rCutSq = rCut * rCut;
  rCutLow = val.ff.cutoffLow;
  rCutLowSq = rCutLow * rCutLow;
  scaling_14 = val.elect.oneFourScale;
  beta = 1 / T_in_K;

  vdwKind = val.ff.VDW_KIND;
  exckind = val.exclude.EXCLUDE_KIND;
  freeEnergy = val.freeEn.enable;

  electrostatic = val.elect.enable;
  ewald = val.elect.ewald;
  wolf = val.elect.wolf;
  coulKind = val.ff.COUL_KIND;
  wolfKind = val.ff.WOLF_KIND;
  tolerance = val.elect.tolerance;
  rswitch = val.ff.rswitch;
  dielectric = val.elect.dielectric;

  if (val.freeEn.enable) {
    sc_alpha = val.freeEn.scaleAlpha;
    sc_sigma = val.freeEn.scaleSigma;
    sc_power = val.freeEn.scalePower;
    sc_coul = val.freeEn.scaleCoulomb;
  } else if (val.neMTMCVal.enable) {
    sc_alpha = val.neMTMCVal.scaleAlpha;
    sc_sigma = val.neMTMCVal.scaleSigma;
    sc_power = val.neMTMCVal.scalePower;
    sc_coul = val.neMTMCVal.scaleCoulomb;
  } else {
    sc_alpha = 0.0;
    sc_sigma = 0.0;
    sc_power = 0;
    sc_coul = false;
  }
  sc_sigma_6 = pow(sc_sigma, 6.0);

  for (uint b = 0; b < BOX_TOTAL; b++) {
    rCutCoulomb[b] = val.elect.cutoffCoulomb[b];
    rCutCoulombSq[b] = rCutCoulomb[b] * rCutCoulomb[b];
    alpha[b] = sqrt(-log(tolerance)) / rCutCoulomb[b];
    alphaSq[b] = alpha[b] * alpha[b];
    recip_rcut[b] = -2.0 * log(tolerance) / rCutCoulomb[b];
    recip_rcut_Sq[b] = recip_rcut[b] * recip_rcut[b];
    if (wolf){
      wolfAlpha[b] = val.elect.wolfAlpha[b];
      wolfFactor1[b] = erfc(wolfAlpha[b]*rCutCoulomb[b])/rCutCoulomb[b];
      wolfFactor2[b] = wolfFactor1[b]/rCutCoulomb[b];
      wolfFactor2[b] += wolfAlpha[b] *  M_2_SQRTPI * 
                        exp(-1.0*wolfAlpha[b]*wolfAlpha[b]*rCutCoulombSq[b])
                        /rCutCoulomb[b];
      wolfFactor3[b] = wolfAlpha[b] *  M_2_SQRTPI;
    }
  }

  vdwGeometricSigma = val.ff.vdwGeometricSigma;
  isMartini = ffKind.isMARTINI;
  exp6 = (vdwKind == val.ff.VDW_EXP6_KIND);

#if ENSEMBLE == GCMC
  isFugacity = val.chemPot.isFugacity;
#endif

  if (vdwKind == val.ff.VDW_STD_KIND)
    particles = new FFParticle(*this);
  else if (vdwKind == val.ff.VDW_EXP6_KIND)
    particles = new FF_EXP6(*this);
  else if (vdwKind == val.ff.VDW_SHIFT_KIND)
    particles = new FF_SHIFT(*this);
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && ffKind.isMARTINI)
    particles = new FF_SWITCH_MARTINI(*this);
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && !ffKind.isMARTINI)
    particles = new FF_SWITCH(*this);
  else {
    std::cout << "Undefined Potential Type detected!\n"
              << "Exiting!\n";
    exit(EXIT_FAILURE);
  }

  if (ffKind.isMARTINI)
    angles = new FFAngleMartini();
  else
    angles = new FFAngles();

  // Define type of interaction to be included. ex. 1-3, 1-4 and more
  if (exckind == val.exclude.EXC_ONETWO_KIND) {
    OneThree = true, OneFour = true, OneN = true;
  } else if (exckind == val.exclude.EXC_ONETHREE_KIND) {
    OneThree = false, OneFour = true, OneN = true;
  } else if (exckind == val.exclude.EXC_ONEFOUR_KIND) {
    OneThree = false, OneFour = false, OneN = true;
  } else {
    std::cout << "Error: Unknown exclude value.\n";
    exit(EXIT_FAILURE);
  }
}

uint Forcefield::GetWolfKind(void) {return wolfKind;}
uint Forcefield::GetCoulKind(void) {return coulKind;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetWolfAlpha(void) {return wolfAlpha;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetWolfFactor1(void) {return wolfFactor1;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetWolfFactor2(void) {return wolfFactor2;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetWolfFactor3(void) {return wolfFactor3;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetRCutCoulomb(void) {return rCutCoulomb;} 
BOX_SIZE_DOUBLE_ARRAY& Forcefield::GetRCutCoulombSq(void) {return rCutCoulombSq;}   
void Forcefield::SetWolfKind(uint wk){
  wolfKind = wk;
}
void Forcefield::SetCoulKind(uint ck){
  coulKind = ck;
}
void Forcefield::SetWolfAlphaAndWolfFactors(double wa, uint b){
    wolfAlpha[b] = wa;
    wolfFactor1[b] = erfc(wolfAlpha[b]*rCutCoulomb[b])/rCutCoulomb[b];
    wolfFactor2[b] = wolfFactor1[b]/rCutCoulomb[b];
    wolfFactor2[b] += wolfAlpha[b] *  M_2_SQRTPI * 
                      exp(-1.0*wolfAlpha[b]*wolfAlpha[b]*rCutCoulombSq[b])
                      /rCutCoulomb[b];
    wolfFactor3[b] = wolfAlpha[b] *  M_2_SQRTPI;
}