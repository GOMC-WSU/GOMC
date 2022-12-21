/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "Wolf.h"

#include <cassert>

#include "BasicTypes.h" //uint
#include "BoxDimensions.h"
#include "CalculateEnergy.h"
#include "Coordinates.h"
#include "EnergyTypes.h"          //Energy structs
#include "EnsemblePreprocessor.h" //Flags
#include "Forcefield.h"
#include "GeomLib.h"
#include "MoleculeKind.h"
#include "NumLib.h"
#include "StaticVals.h" //For init
#include "System.h"     //For init
#include "TrialMol.h"
#ifdef GOMC_CUDA
#include "CalculateEwaldCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
#include "GOMCEventsProfile.h"

//
//
//    Energy Calculation functions for Wolf approximation method
//    Calculating self and correction terms.
//
//    Developed by G. Schwing and Mohammad S. Barhaghi
//
//

using namespace geom;

Wolf::Wolf(StaticVals &stat, System &sys)
    : ff(stat.forcefield), mols(stat.mol), currentCoords(sys.coordinates),
#ifdef VARIABLE_PARTICLE_NUMBER
      molLookup(sys.molLookup),  
#else
      molLookup(stat.molLookup),
#endif
      currentAxes(sys.boxDimRef), currentCOM(sys.com), sysPotRef(sys.potential),
      lambdaRef(sys.lambdaRef),
      wolfKind(stat.forcefield.wolfKind),
      coulKind(stat.forcefield.coulKind),
      wolfAlpha(stat.forcefield.wolfAlpha),
      wolfFactor1(stat.forcefield.wolfFactor1),      
      wolfFactor2(stat.forcefield.wolfFactor2),
      wolfFactor3(stat.forcefield.wolfFactor3),
      rCutCoulomb(stat.forcefield.rCutCoulomb),
      rCutCoulombSq(stat.forcefield.rCutCoulombSq)  {
  ewald = false;
  electrostatic = false;
  alpha = 0.0;
  recip_rcut = 0.0;
  recip_rcut_Sq = 0.0;
  multiParticleEnabled = stat.multiParticleEnabled;
  oneThree = stat.forcefield.OneThree;
  oneFour = stat.forcefield.OneFour;
  scaling_14 = stat.forcefield.scaling_14;
}

void Wolf::Init() {
  for (uint m = 0; m < mols.count; ++m) {
    const MoleculeKind &molKind = mols.GetKind(m);
    for (uint a = 0; a < molKind.NumAtoms(); ++a) {
      particleKind.push_back(molKind.AtomKind(a));
      particleMol.push_back(m);
      particleCharge.push_back(molKind.AtomCharge(a));
      if (std::abs(molKind.AtomCharge(a)) < 0.000000001) {
        particleHasNoCharge.push_back(true);
      } else {
        particleHasNoCharge.push_back(false);
      }
    }
  }

  // initialize starting index and length index of each molecule
  startMol.resize(currentCoords.Count());
  lengthMol.resize(currentCoords.Count());

  for (int atom = 0; atom < (int)currentCoords.Count(); atom++) {
    startMol[atom] = mols.MolStart(particleMol[atom]);
    lengthMol[atom] = mols.MolLength(particleMol[atom]);
  }
}


// calculate correction term for a molecule, with system lambda
double Wolf::MolCorrection(uint molIndex, uint box) const {
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_MOL);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0, undampenedCorr = 0.0;
  XYZ virComponents;

  MoleculeKind &thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      correction = RahbariCorrection(thisKind,
                              box,
                              virComponents,
                              atomSize,
                              start);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      correction = Waibel2018Correction(thisKind,
                              box,
                              virComponents,
                              atomSize,
                              start); 
      break; 
    case WOLF_WAIBEL2019_KIND :
      correction = Waibel2019Correction(thisKind,
                        box,
                        virComponents,
                        atomSize,
                        start);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_MOL);
  return correction * lambdaCoef * lambdaCoef;
}

// It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void Wolf::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                             const std::vector<double> &lambda_Coul,
                             const uint iState, const uint molIndex,
                             const uint box) const {
  uint atomSize = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  uint lambdaSize = lambda_Coul.size();
  double coefDiff, distSq, dist, correction = 0.0;
  XYZ virComponents;
  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];
  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      correction = RahbariCorrection(thisKind,
                              box,
                              virComponents,
                              atomSize,
                              start);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      correction = Waibel2018Correction(thisKind,
                              box,
                              virComponents,
                              atomSize,
                              start); 
      break; 
    case WOLF_WAIBEL2019_KIND :
      correction = Waibel2019Correction(thisKind,
                        box,
                        virComponents,
                        atomSize,
                        start);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);  
  }
  // Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].correction += coefDiff * correction;
  }
  // Calculate du/dl of correction for current state, for linear scaling
  dUdL_Coul.correction += correction;
}

// calculate self term for a box, using system lambda
double Wolf::BoxSelf(uint box) const {
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::SELF_BOX);
  double self = 0.0;
  double molSelfEnergy;
  uint i, j, length, molNum;
  double lambdaCoef = 1.0;

  for (i = 0; i < mols.GetKindsCount(); i++) {
    MoleculeKind const &thisKind = mols.kinds[i];
    length = thisKind.NumAtoms();
    molNum = molLookup.NumKindInBox(i, box);
    molSelfEnergy = 0.0;
    if (lambdaRef.KindIsFractional(i, box)) {
      // If a molecule is fractional, we subtract the fractional molecule and
      // add it later
      --molNum;
      // returns lambda and not sqrt(lambda)
      lambdaCoef = lambdaRef.GetLambdaCoulomb(i, box);
    }

    for (j = 0; j < length; j++) {
      molSelfEnergy += (thisKind.AtomCharge(j) * thisKind.AtomCharge(j));
    }
    self += (molSelfEnergy * molNum);
    if (lambdaRef.KindIsFractional(i, box)) {
      // Add the fractional molecule part
      self += (molSelfEnergy * lambdaCoef);
    }
  }

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      // we eliminate the alpha/root(pi) using Wolf,mod from WAIBEL2018 et al
      self *= -0.5 * wolfFactor1[box];
      break; 
    case WOLF_WAIBEL2019_KIND :
      self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_BOX);
  return self * num::qqFact;
}


// calculate correction term for a molecule with lambda = 1
// It's called when the molecule configuration changes, moleculeTransfer, MEMC
// It never been called in Free Energy calculation, because we are in
// NVT and NPT ensemble
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double correction = 0.0;

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      correction = RahbariCorrection(trialMol);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      correction = Waibel2018Correction(trialMol);
      break; 
    case WOLF_WAIBEL2019_KIND :
      correction = Waibel2019Correction(trialMol);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);  
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return correction;
}

// calculate correction term for a molecule with system lambda
// It's called when the molecule configuration changes, regrowth, crankshaft,
// IntraSwap, IntraMEMC ...
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol,
                             const uint molIndex) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double lambdaCoef = GetLambdaCoef(molIndex, box);
  double correction = 0.0;

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      correction = RahbariCorrection(trialMol);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      correction = Waibel2018Correction(trialMol);
      break; 
    case WOLF_WAIBEL2019_KIND :
      correction = Waibel2019Correction(trialMol);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);  
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return correction * lambdaCoef * lambdaCoef;
}

// It's called if we transfer one molecule from one box to another
// No need to scale the charge with lambda, since this function is not being
// called from free energy or NeMTMC
double Wolf::SwapSelf(const cbmc::TrialMol &trialMol) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::SELF_SWAP);
  MoleculeKind const &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  double en_self = 0.0;

  for (uint i = 0; i < atomSize; i++) {
    en_self += (thisKind.AtomCharge(i) * thisKind.AtomCharge(i));
  }

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      en_self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      // we eliminate the alpha/root(pi) using Wolf,mod from WAIBEL2018 et al
      en_self *= -0.5 * wolfFactor1[box];
      break; 
    case WOLF_WAIBEL2019_KIND :
      en_self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_SWAP);
  return en_self * num::qqFact;
}

// It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                       const std::vector<double> &lambda_Coul,
                       const uint iState, const uint molIndex,
                       const uint box) const {
  uint atomSize = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  uint lambdaSize = lambda_Coul.size();
  double coefDiff, en_self = 0.0;
  // Calculate the self energy with lambda = 1
  for (uint i = 0; i < atomSize; i++) {
    en_self += (particleCharge[i + start] * particleCharge[i + start]);
  }

  switch(wolfKind) {
    case WOLF_RAHBARI_KIND  :
      en_self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break; 
    case WOLF_WAIBEL2018_KIND  :
      // we eliminate the alpha/root(pi) using Wolf,mod from WAIBEL2018 et al
      en_self *= -0.5 * wolfFactor1[box];
      break; 
    case WOLF_WAIBEL2019_KIND :
      en_self *= -0.5 * ((wolfAlpha[box] * M_2_SQRTPI) + wolfFactor1[box]);
      break;
    // you can have any number of case statements.
    default : //Optional
      exit(1);
  }
  en_self *= num::qqFact;
  
  // Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].self += coefDiff * en_self;
  }
  // Calculate du/dl of self for current state, for linear scaling
  dUdL_Coul.self += en_self;
}


double Wolf::GetLambdaCoef(uint molA, uint box) const {
  double lambda = lambdaRef.GetLambdaCoulomb(molA, box);
  // Each charge gets sq root of it.
  return sqrt(lambda);
}

double Wolf::RahbariCorrection(MoleculeKind &thisKind,
                                    uint box,
                                    XYZ virComponents,
                                    uint atomSize,
                                    uint start) const{
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;
  for (uint i = 0; i < atomSize; i++) {
    if (particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords, start + i,
                        start + j, box);
      dist = sqrt(distSq);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq[box]){
        // All methods have this constant term.
        dampenedCorr -= erf(wolfAlpha[box] * dist)/dist;   
        dampenedCorr -= wolfFactor1[box];
        // Rahbari doesn't include this term
        //if(coulKind){
        //  double distDiff = dist-rCutCoulomb[box];
        //  dampenedCorr += wolfFactor2[box]*distDiff;
        //} 
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
        //correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
        //              erf(ff.alpha[box] * dist) / dist);
      }
    }
  }
  return correction * num::qqFact;
}


double Wolf::Waibel2018Correction(MoleculeKind &thisKind,
                                    uint box,
                                    XYZ virComponents,
                                    uint atomSize,
                                    uint start) const{
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;
  for (uint i = 0; i < atomSize; i++) {
    if (particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      correction -= thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * wolfFactor1[box];
    }
    if(oneThree){
      //loop over all 1-3 partners of the particle
      const uint* partner = thisKind.sortedNB_1_3.Begin(i);
      const uint* end = thisKind.sortedNB_1_3.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
          dampenedCorr *= scaling_14;
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }
    }
    if(oneFour){
      //loop over all 1-4 partners of the particle
      const uint* partner = thisKind.sortedNB_1_4.Begin(i);
      const uint* end = thisKind.sortedNB_1_4.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
          dampenedCorr *= scaling_14;
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }
    }
    //loop over all 1-N partners of the particle
    const uint* partner = thisKind.sortedNB.Begin(i);
    const uint* end = thisKind.sortedNB.End(i);
    while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
      }
      ++partner;
    }      
  }
  return correction * num::qqFact;
}

double Wolf::Waibel2019Correction(MoleculeKind &thisKind,
                                    uint box,
                                    XYZ virComponents,
                                    uint atomSize,
                                    uint start) const{
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;
  for (uint i = 0; i < atomSize; i++) {
    if (particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords, start + i,
                        start + j, box);
      dist = sqrt(distSq);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq[box]){
        // All methods have this constant term.
        dampenedCorr -= erf(wolfAlpha[box] * dist)/dist;   
        dampenedCorr -= wolfFactor1[box];
        if(coulKind){
          double distDiff = dist-rCutCoulomb[box];
          dampenedCorr += wolfFactor2[box]*distDiff;
        } 
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
        //correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
        //              erf(ff.alpha[box] * dist) / dist);
      }
    }
  }
  return correction * num::qqFact;
}


double Wolf::RahbariCorrection(const cbmc::TrialMol &trialMol) const{
  uint box = trialMol.GetBox();
  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  XYZ virComponents;
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;
  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq[box]){
        // All methods have this constant term.
        dampenedCorr -= wolfFactor1[box];
        dist = sqrt(distSq);
        dampenedCorr += -1.0*erf(wolfAlpha[box] * dist)/dist;   
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  return correction * num::qqFact;
}


double Wolf::Waibel2018Correction(const cbmc::TrialMol &trialMol) const{
  uint box = trialMol.GetBox();
  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  XYZ virComponents;
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;
  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      correction -= thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * wolfFactor1[box];
    }
    if(oneThree){
      //loop over all 1-3 partners of the particle
      const uint* partner = thisKind.sortedNB_1_3.Begin(i);
      const uint* end = thisKind.sortedNB_1_3.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                        i, *partner, box); 
        if(distSq < rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
          dampenedCorr *= scaling_14;
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }
    }
    if(oneFour){
      //loop over all 1-4 partners of the particle
      const uint* partner = thisKind.sortedNB_1_4.Begin(i);
      const uint* end = thisKind.sortedNB_1_4.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                        i, *partner, box); 
        if(distSq < rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
          dampenedCorr *= scaling_14;
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }
    }
    //loop over all 1-N partners of the particle
    const uint* partner = thisKind.sortedNB.Begin(i);
    const uint* end = thisKind.sortedNB.End(i);
    while (partner != end) {
      // Need to check for cutoff for all kinds
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                      i, *partner, box); 
      if(distSq < rCutCoulombSq[box]){
        dist = sqrt(distSq);
        dampenedCorr = -1.0*erf(wolfAlpha[box] * dist)/dist;  
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
      }
      ++partner;
    }      
  }
  return correction * num::qqFact;
}

double Wolf::Waibel2019Correction(const cbmc::TrialMol &trialMol) const{
  uint box = trialMol.GetBox();
  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  XYZ virComponents;
  double correction = 0.0, dampenedCorr = 0.0;
  double dist, distSq;

  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq[box]){
        // All methods have this constant term.
        dampenedCorr -= wolfFactor1[box];
        dist = sqrt(distSq);
        dampenedCorr += -1.0*erf(wolfAlpha[box] * dist)/dist;   
        if(coulKind){
          double distDiff = dist-rCutCoulomb[box];
          dampenedCorr += wolfFactor2[box]*distDiff;
        }
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  return correction * num::qqFact;
}