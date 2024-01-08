/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef EWALD_CACHED_H
#define EWALD_CACHED_H

#include "Ewald.h"

class EwaldCached : public Ewald {
public:
  EwaldCached(StaticVals &stat, System &sys);
  virtual ~EwaldCached();

  virtual void Init();

  virtual void AllocMem();

  // compute reciprocal term for a box with a new volume
  virtual void BoxReciprocalSetup(uint box, XYZArray const &molCoords);

  // compute reciprocal term for a box when not testing a volume change
  virtual void BoxReciprocalSums(uint box, XYZArray const &molCoords);

  // calculate reciprocal energy term for a box
  virtual double BoxReciprocal(uint box, bool isNewVolume) const;

  // calculate reciprocal term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                               const uint box);

  // calculate reciprocal term for lambdaNew and Old with same coordinates
  virtual double ChangeLambdaRecip(XYZArray const &molCoords,
                                   const double lambdaOld,
                                   const double lambdaNew, const uint molIndex,
                                   const uint box);

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
                        const std::vector<uint> &molIndexOld);

  // It's called in free energy calculation to calculate the change in
  // reciprocal energy in all lambda states
  virtual void ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                           const std::vector<double> &lambda_Coul,
                           const uint iState, const uint molIndex,
                           const uint box) const;
  // restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  // update sinMol and cosMol
  virtual void exgMolCache();

  // backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache();

private:
  double *cosMolRestore; // cos()*charge
  double *sinMolRestore; // sin()*charge
  double **cosMolRef;
  double **sinMolRef;
  double **cosMolBoxRecip;
  double **sinMolBoxRecip;
#if ENSEMBLE == GEMC
  const uint GEMC_KIND;
#endif
};

#endif /*EWALD_CACHED_H*/
