/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef NOEWALD_H
#define NOEWALD_H

#include "Ewald.h"

class NoEwald : public Ewald {
  // friend class CalculateEnergy;
public:
  NoEwald(StaticVals &stat, System &sys);

  virtual void Init();

  virtual void AllocMem();

  // initialize term used for ewald calculation
  virtual void RecipInit(uint box, BoxDimensions const &boxAxes);

  // calculate self term for a box
  virtual double BoxSelf(uint box) const;

  // compute reciprocal term for a box with a new volume
  virtual void BoxReciprocalSetup(uint box, XYZArray const &molCoords);

  // compute reciprocal term for a box when not testing a volume change
  virtual void BoxReciprocalSums(uint box, XYZArray const &molCoords);

  // calculate reciprocal energy term for a box
  virtual double BoxReciprocal(uint box, bool isNewVolume) const;

  // calculate reciprocal force term for a box with molCoords
  virtual void BoxForceReciprocal(XYZArray const &molCoords,
                                  XYZArray &atomForceRec, XYZArray &molForceRec,
                                  uint box);

  // calculate reciprocal force term for a box
  virtual Virial VirialReciprocal(Virial &virial, uint box) const;

  // calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box) const;

  // calculate reciprocal term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                               const uint box);

  // calculate reciprocal term for lambdaNew and Old with same coordinates
  virtual double ChangeLambdaRecip(XYZArray const &molCoords,
                                   const double lambdaOld,
                                   const double lambdaNew, const uint molIndex,
                                   const uint box);

  // calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol &trialMo) const;

  // calculate correction term after swap move with lambda = 1
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol) const;

  // calculate correction term after swap move, with system lambda
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol,
                                const uint molIndex) const;

  // calculate reciprocal term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex);

  // calculate reciprocal term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box,
                                 const int molIndex);

  // calculate reciprocal term for inserting some molecules (kindA) in
  // destination box and removing a molecule (kindB) from destination box
  virtual double
  MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                        const std::vector<cbmc::TrialMol> &oldMol,
                        const std::vector<uint> &molIndexNew,
                        const std::vector<uint> &molIndexold, bool first_call);

  // back up reciprocal value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box);

  // It's called in free energy calculation to calculate the change in
  // self energy in all lambda states
  virtual void ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const;

  // It's called in free energy calculation to calculate the change in
  // correction energy in all lambda states
  virtual void ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                                const std::vector<double> &lambda_Coul,
                                const uint iState, const uint molIndex,
                                const uint box) const;

  // It's called in free energy calculation to calculate the change in
  // reciprocal energy in all lambda states
  virtual void ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                           const std::vector<double> &lambda_Coul,
                           const uint iState, const uint molIndex,
                           const uint box) const;

  // update reciprocal values
  virtual void UpdateRecip(uint box);

  // copy reciprocal values from ref to new
  virtual void CopyRecip(uint box);

  // update the kx, ky, kz, hsqr and prefact
  virtual void UpdateRecipVec(uint box);

  // restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  // update sinMol and cosMol
  virtual void exgMolCache();

  // backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache();

  virtual void UpdateVectorsAndRecipTerms(bool output);
};

#endif /*NOEWALD_H*/
