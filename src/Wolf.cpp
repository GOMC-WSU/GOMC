/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "Wolf.h"

Wolf::Wolf(StaticVals &stat, System &sys) : Ewald(stat, sys) {}

void Wolf::Init() {}

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
double Wolf::BoxSelf(uint box) const { return 0.0; }

// calculate correction term for a molecule
double Wolf::MolCorrection(uint molIndex, uint box) const { return 0.0; }

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

// calculate self term after swap move
double Wolf::SwapSelf(const cbmc::TrialMol &trialMol) const { return 0.0; }

// calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol) const {
  return 0.0;
}
// calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol &trialMol,
                               const uint molIndex) const {
  return 0.0;
}

// It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box) const {
  return;
}

// It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void Wolf::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                               const std::vector<double> &lambda_Coul,
                               const uint iState, const uint molIndex,
                               const uint box) const {
  return;
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
