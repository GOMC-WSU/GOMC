/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "NoEwald.h"

NoEwald::NoEwald(StaticVals &stat, System &sys) : Ewald(stat, sys) {}

void NoEwald::Init() {}

void NoEwald::AllocMem() { return; }

void NoEwald::RecipInit(uint box, BoxDimensions const &boxAxes) { return; }

// compute reciprocal term for a box with a new volume
void NoEwald::BoxReciprocalSetup(uint box, XYZArray const &molCoords) {
  return;
}

// compute reciprocal term for a box when not testing a volume change
void NoEwald::BoxReciprocalSums(uint box, XYZArray const &molCoords) { return; }

// calculate reciprocal term for a box
double NoEwald::BoxReciprocal(uint box, bool isNewVolume) const { return 0.0; }

// calculate reciprocal force term for a box with molCoords
void NoEwald::BoxForceReciprocal(XYZArray const &molCoords,
                                 XYZArray &atomForceRec, XYZArray &molForceRec,
                                 uint box) {
  return;
}

// calculate reciprocal force term for a box
Virial NoEwald::VirialReciprocal(Virial &virial, uint box) const {
  return virial;
}

// calculate reciprocal term for displacement and rotation move
double NoEwald::MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                              const uint box) {
  return 0.0;
}

// calculate reciprocal term for lambdaNew and Old with same coordinates
double NoEwald::ChangeLambdaRecip(XYZArray const &molCoords,
                                  const double lambdaOld,
                                  const double lambdaNew, const uint molIndex,
                                  const uint box) {
  return 0.0;
}

// calculate self term for a box
double NoEwald::BoxSelf(uint box) const { return 0.0; }

// calculate correction term for a molecule
double NoEwald::MolCorrection(uint molIndex, uint box) const { return 0.0; }

// calculate reciprocal term in destination box for swap move
double NoEwald::SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                              const int molIndex) {
  return 0.0;
}

// calculate reciprocal term in source box for swap move
double NoEwald::SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box,
                                const int molIndex) {
  return 0.0;
}

// calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double NoEwald::MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                                      const std::vector<cbmc::TrialMol> &oldMol,
                                      const std::vector<uint> &molIndexNew,
                                      const std::vector<uint> &molIndexold,
                                      bool first_call) {
  return 0.0;
}

// calculate self term after swap move
double NoEwald::SwapSelf(const cbmc::TrialMol &trialMol) const { return 0.0; }

// calculate correction term after swap move
double NoEwald::SwapCorrection(const cbmc::TrialMol &trialMol) const {
  return 0.0;
}
// calculate correction term after swap move
double NoEwald::SwapCorrection(const cbmc::TrialMol &trialMol,
                               const uint molIndex) const {
  return 0.0;
}

// It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void NoEwald::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box) const {
  return;
}

// It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void NoEwald::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                               const std::vector<double> &lambda_Coul,
                               const uint iState, const uint molIndex,
                               const uint box) const {
  return;
}

// It's called in free energy calculation to calculate the change in
// reciprocal energy in all lambda states
void NoEwald::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const {
  return;
}

// back up reciprocal value to Ref (will be called during initialization)
void NoEwald::SetRecipRef(uint box) { return; }

// update reciprocal values
void NoEwald::UpdateRecip(uint box) { return; }

// copy reciprocal values from ref to new
void NoEwald::CopyRecip(uint box) { return; }

// update the kx, ky, kz, hsqr and prefact
void NoEwald::UpdateRecipVec(uint box) { return; }

// restore cosMol and sinMol
void NoEwald::RestoreMol(int molIndex) { return; }

// restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::exgMolCache() { return; }

// backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::backupMolCache() { return; }

void NoEwald::UpdateVectorsAndRecipTerms(bool output) { return; }
