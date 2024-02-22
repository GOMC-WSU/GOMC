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
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocal part of ewald
//
//    Developed by Y. Li and Mohammad S. Barhaghi
//
//

using namespace geom;

Wolf::Wolf(StaticVals &stat, System &sys) : Ewald(stat, sys) {}
Wolf::~Wolf() {}

//void Wolf::Init() {}
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

  //AllocMem();
  // initialize K vectors and reciprocal terms
  //UpdateVectorsAndRecipTerms(true);
}

void Wolf::AllocMem() { return; }

void Wolf::RecipInit(uint box, BoxDimensions const &boxAxes) { return; }

// compute reciprocal term for a box with a new volume
void Wolf::BoxReciprocalSetup(uint box, XYZArray const &molCoords) {
  return;
}

// compute reciprocal term for a box when not testing a volume change
void Wolf::BoxReciprocalSums(uint box, XYZArray const &molCoords) { return; }

// calculate reciprocal term for a box
double Wolf::BoxReciprocal(uint box, bool isNewVolume) const { return 0.0; }

// calculate reciprocal force term for a box with molCoords
void Wolf::BoxForceReciprocal(XYZArray const &molCoords,
                                 XYZArray &atomForceRec, XYZArray &molForceRec,
                                 uint box) {
  return;
}

// calculate reciprocal force term for a box
Virial Wolf::VirialReciprocal(Virial &virial, uint box) const {
  return virial;
}

// calculate reciprocal term for displacement and rotation move
double Wolf::MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                              const uint box) {
  return 0.0;
}

// calculate reciprocal term for lambdaNew and Old with same coordinates
double Wolf::ChangeLambdaRecip(XYZArray const &molCoords,
                                  const double lambdaOld,
                                  const double lambdaNew, const uint molIndex,
                                  const uint box) {
  return 0.0;
}

// calculate self term for a box
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
  if (ff.simpleself)
    self *= -0.5 * (ff.wolf_factor_1[box]) * num::qqFact;
  else
    self *= -0.5 * ((ff.wolf_alpha[box] * M_2_SQRTPI) + ff.wolf_factor_1[box]) * num::qqFact;
  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_BOX);
  return self;
}

// calculate correction term for a molecule, with system lambda
double Wolf::MolCorrection(uint molIndex, uint box) const {
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_MOL);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;

  MoleculeKind &thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);
  
  if(ff.simpleself){
    correction = SimpleSelfCorrection(molIndex, box);
  } else {
    for (uint i = 0; i < atomSize; i++) {
      if (particleHasNoCharge[start + i]) {
        continue;
      }
      for (uint j = i + 1; j < atomSize; j++) {
        currentAxes.InRcut(distSq, virComponents, currentCoords, start + i,
                          start + j, box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          //correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
          //               erf(ff.alpha[box] * dist) / dist);
          correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            correction += (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * ff.wolf_factor_2[box]*distDiff);
          }
        }
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_MOL);
  return -1.0 * num::qqFact * correction * lambdaCoef * lambdaCoef;
}


// calculate reciprocal term in destination box for swap move
double Wolf::SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                              const int molIndex) {
  return 0.0;
}

// calculate reciprocal term in source box for swap move
double Wolf::SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box,
                                const int molIndex) {
  return 0.0;
}

// calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double Wolf::MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                                      const std::vector<cbmc::TrialMol> &oldMol,
                                      const std::vector<uint> &molIndexNew,
                                      const std::vector<uint> &molIndexold,
                                      bool first_call) {
  return 0.0;
}

// calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr= 0.0;
  XYZ virComponents;
  const MoleculeKind &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();

  if(ff.simpleself){
    correction = SimpleSelfCorrection(trialMol);
  } else {
    for (uint i = 0; i < atomSize; i++) {
      for (uint j = i + 1; j < atomSize; j++) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, j,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          //correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
          //               erf(ff.alpha[box] * dist) / dist);
          correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * ff.wolf_factor_2[box]*distDiff);
          }
        }
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return num::qqFact * correction;
}

// calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol,
                               const uint molIndex) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;
  const MoleculeKind &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);

  if(ff.simpleself){
    correction = SimpleSelfCorrection(trialMol,molIndex);
  } else {
    for (uint i = 0; i < atomSize; i++) {
      if (particleHasNoCharge[start + i]) {
        continue;
      }
      for (uint j = i + 1; j < atomSize; j++) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, j,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          //correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
          //               erf(ff.alpha[box] * dist) / dist);
          correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            correction -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * ff.wolf_factor_2[box]*distDiff);
          }
        }
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return num::qqFact * correction * lambdaCoef * lambdaCoef;
}


// calculate self term after swap move
double Wolf::SwapSelf(const cbmc::TrialMol &trialMol) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::SELF_SWAP);
  MoleculeKind const &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  double en_self = 0.0;

  for (uint i = 0; i < atomSize; i++) {
    en_self -= (thisKind.AtomCharge(i) * thisKind.AtomCharge(i));
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_SWAP);
  // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
  //return (en_self * ff.alpha[box] * num::qqFact * M_2_SQRTPI * 0.5);
  if (ff.simpleself)
    return (en_self *  num::qqFact * (ff.wolf_factor_1[box]));
  else
    return (en_self *  num::qqFact * ((ff.wolf_alpha[box]*M_2_SQRTPI * 0.5) + ff.wolf_factor_1[box]));
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
  // M_2_SQRTPI is 2/sqrt(PI), so need to multiply by 0.5 to get sqrt(PI)
  //en_self *= -1.0 * ff.alpha[box] * num::qqFact * M_2_SQRTPI * 0.5;
  if (ff.simpleself)
    en_self *= -1.0 * num::qqFact * (ff.wolf_factor_1[box]);
  else
    en_self *= -1.0 * num::qqFact * ((ff.wolf_alpha[box]*M_2_SQRTPI * 0.5) + ff.wolf_factor_1[box]);

  // Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].self += coefDiff * en_self;
  }
  // Calculate du/dl of self for current state, for linear scaling
  dUdL_Coul.self += en_self;
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
  MoleculeKind &thisKind = mols.kinds[mols.kIndex[molIndex]];
  double coefDiff, distSq, dist, correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;

  // Calculate the correction energy with lambda = 1
  if (ff.simpleself){
    correction = SimpleSelfCorrection(molIndex,box);
  } else {
    for (uint i = 0; i < atomSize; i++) {
      if (particleHasNoCharge[start + i]) {
        continue;
      }

      for (uint j = i + 1; j < atomSize; j++) {
        distSq = 0.0;
        currentAxes.InRcut(distSq, virComponents, currentCoords, start + i,
                          start + j, box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          //correction += (particleCharge[i + start] * particleCharge[j + start] *
          //               erf(ff.alpha[box] * dist) / dist);
          correction += (particleCharge[i + start] * particleCharge[j + start] *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            correction += (particleCharge[i + start] * particleCharge[j + start] * ff.wolf_factor_2[box]*distDiff);
          }
        }
      }
    }
  }
  correction *= -1.0 * num::qqFact;
  // Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].correction += coefDiff * correction;
  }
  // Calculate du/dl of correction for current state, for linear scaling
  dUdL_Coul.correction += correction;
}

// It's called in free energy calculation to calculate the change in
// reciprocal energy in all lambda states
void Wolf::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const {
  return;
}

// back up reciprocal value to Ref (will be called during initialization)
void Wolf::SetRecipRef(uint box) { return; }

// update reciprocal values
void Wolf::UpdateRecip(uint box) { return; }

// copy reciprocal values from ref to new
void Wolf::CopyRecip(uint box) { return; }

// update the kx, ky, kz, hsqr and prefact
void Wolf::UpdateRecipVec(uint box) { return; }

// restore cosMol and sinMol
void Wolf::RestoreMol(int molIndex) { return; }

// restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::exgMolCache() { return; }

// backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::backupMolCache() { return; }

void Wolf::UpdateVectorsAndRecipTerms(bool output) { return; }

double Wolf::SimpleSelfCorrection(const uint molIndex,
                                  const uint box) const {
  uint atomSize = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  MoleculeKind &thisKind = mols.kinds[mols.kIndex[molIndex]];
  double coefDiff, distSq, dist, correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;
  for (uint i = 0; i < atomSize; i++) {
    if (particleHasNoCharge[start + i]) {
      continue;
    }
    if(ff.OneThree){
      //loop over all 1-3 partners of the particle
      const uint* partner = thisKind.sortedNB_1_3.Begin(i);
      const uint* end = thisKind.sortedNB_1_3.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < ff.rCutCoulombSq[box]){
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction += (ff.scaling_14*dampenedCorr);
        }
        ++partner;
      }
    }
    if(ff.OneFour){
      //loop over all 1-4 partners of the particle
      const uint* partner = thisKind.sortedNB_1_4.Begin(i);
      const uint* end = thisKind.sortedNB_1_4.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < ff.rCutCoulombSq[box]){
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction += (ff.scaling_14*dampenedCorr);
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
      if(distSq < ff.rCutCoulombSq[box]){
        dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                      ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
        if(ff.intramoleculardsf){
          double distDiff = dist-ff.rCutCoulomb[box];
          dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
        }
        // Dont scale 1-N
        correction += (dampenedCorr);
      }
      ++partner;
    }
  }
  return correction;
}


double Wolf::SimpleSelfCorrection(const cbmc::TrialMol &trialMol,
                               const uint molIndex) const {
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0;
  XYZ virComponents;
  const MoleculeKind &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);
  for (uint i = 0; i < atomSize; i++) {
    if(ff.OneThree){
      //loop over all 1-3 partners of the particle
      const uint* partner = thisKind.sortedNB_1_3.Begin(i);
      const uint* end = thisKind.sortedNB_1_3.End(i);
      while (partner != end) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (ff.scaling_14*dampenedCorr);
        }
        ++partner;
      }
    }
    if(ff.OneFour){
      //loop over all 1-4 partners of the particle
      const uint* partner = thisKind.sortedNB_1_4.Begin(i);
      const uint* end = thisKind.sortedNB_1_4.End(i);
      while (partner != end) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (ff.scaling_14*dampenedCorr);
        }
        ++partner;
      }
    }
    //loop over all 1-4 partners of the particle
    const uint* partner = thisKind.sortedNB.Begin(i);
    const uint* end = thisKind.sortedNB.End(i);
    while (partner != end) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                        box);
      if(distSq < ff.rCutCoulombSq[box]){
        dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (dampenedCorr);
      }
      ++partner;
    }
  }
  return correction;
}


double Wolf::SimpleSelfCorrection(const cbmc::TrialMol &trialMol) const {
  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr= 0.0;
  XYZ virComponents;
  const MoleculeKind &thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  for (uint i = 0; i < atomSize; i++) {
    if(ff.OneThree){
      //loop over all 1-3 partners of the particle
      const uint* partner = thisKind.sortedNB_1_3.Begin(i);
      const uint* end = thisKind.sortedNB_1_3.End(i);
      while (partner != end) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (ff.scaling_14*dampenedCorr);
        }
        ++partner;
      }
    }
    if(ff.OneFour){
      //loop over all 1-4 partners of the particle
      const uint* partner = thisKind.sortedNB_1_4.Begin(i);
      const uint* end = thisKind.sortedNB_1_4.End(i);
      while (partner != end) {
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                          box);
        if(distSq < ff.rCutCoulombSq[box]){
          dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (ff.scaling_14*dampenedCorr);
        }
        ++partner;
      }
    }
    //loop over all 1-4 partners of the particle
    const uint* partner = thisKind.sortedNB.Begin(i);
    const uint* end = thisKind.sortedNB.End(i);
    while (partner != end) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(), i, *partner,
                        box);
      if(distSq < ff.rCutCoulombSq[box]){
        dist = sqrt(distSq);
          dampenedCorr = (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) *
                        ((erf(ff.wolf_alpha[box] * dist) / dist) + ff.wolf_factor_1[box]));
          if(ff.intramoleculardsf){
            double distDiff = dist-ff.rCutCoulomb[box];
            dampenedCorr += (thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * ff.wolf_factor_2[box]*distDiff);
          }
          correction -= (dampenedCorr);
      }
      ++partner;
    }
  }

  return correction;
}