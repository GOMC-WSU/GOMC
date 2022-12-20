/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef ELECTROSTATIC_BASE_H
#define ELECTROSTATIC_BASE_H

#include "EnsemblePreprocessor.h"
#include "TrialMol.h"

// An abstract class
class ElectrostaticBase
{
	// Data members of class
public:
	// Pure Virtual Function
	//virtual void show() = 0;

  virtual void Init() = 0;

  virtual void AllocMem() = 0;

  // initialize term used for ewald calculation
  virtual void RecipInit(uint box, BoxDimensions const &boxAxes) = 0;

  // initialize wave vector for orthogonal box
  virtual void RecipInitOrth(uint box, BoxDimensions const &boxAxes) = 0;

  // initialize wave vector for non-orthogonal box
  virtual void RecipInitNonOrth(uint box, BoxDimensions const &boxAxes) = 0;

  // Get initial estimate of memory required
  virtual void RecipCountInit(uint box, BoxDimensions const &boxAxes) = 0;

  // compute reciprocal term for a box with a new volume
  virtual void BoxReciprocalSetup(uint box, XYZArray const &molCoords) = 0;

  // compute reciprocal term for a box when not testing a volume change
  virtual void BoxReciprocalSums(uint box, XYZArray const &molCoords) = 0;

  // calculate reciprocal energy term for a box
  virtual double BoxReciprocal(uint box, bool isNewVolume) const = 0;

  // calculate self term for a box
  virtual double BoxSelf(uint box) const = 0;

  // calculate reciprocal force term for a box
  virtual Virial VirialReciprocal(Virial &virial, uint box) const = 0;

  // calculate reciprocal term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                               const uint box) = 0;

  // calculate reciprocal term for lambdaNew and Old with same coordinates
  virtual double ChangeLambdaRecip(XYZArray const &molCoords,
                                   const double lambdaOld,
                                   const double lambdaNew, const uint molIndex,
                                   const uint box) = 0;

  // calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box) const = 0;

  // calculate reciprocal term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex) = 0;

  // calculate reciprocal term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box,
                                 const int molIndex) = 0;

  // calculate reciprocal term for inserting some molecules (kindA) in
  // destination box and removing a molecule (kindB) from destination box
  virtual double
  MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                        const std::vector<cbmc::TrialMol> &oldMol,
                        const std::vector<uint> &molIndexNew,
                        const std::vector<uint> &molIndexOld, bool first_call) = 0;

  // calculate correction term after swap move, with lambda = 1
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol) const = 0;

  // calculate correction term after swap move, with system lambda
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol,
                                const uint molIndex) const = 0;

  // back up reciprocal value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box) = 0;

  // update reciprocal values
  virtual void UpdateRecip(uint box) = 0;

  // copy reciprocal values from ref to new
  virtual void CopyRecip(uint box) = 0;

  // update kx, ky, kz, hsqr and prefact
  virtual void UpdateRecipVec(uint box) = 0;

  // calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol &trialMol) const = 0;

  // restore cosMol and sinMol
  virtual void RestoreMol(int molIndex) = 0;

  // Find the largest kvector
  virtual uint findLargeImage() = 0;

  // update sinMol and cosMol
  virtual void exgMolCache() = 0;

  // backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache() = 0;

  /// This function performs three actions:
  /// 1. Initialize k vectors
  /// 2. Run BoxReciprocalSetup to calculate sumRnew and sumInew vectors
  /// 3. Then copy all new vectors to ref vectors
  /// @param bool: whether to print the vector size
  virtual void UpdateVectorsAndRecipTerms(bool output) = 0;

  // calculate reciprocal force term for a box with molCoords
  virtual void BoxForceReciprocal(XYZArray const &molCoords,
                                  XYZArray &atomForceRec, XYZArray &molForceRec,
                                  uint box) = 0;

  double GetLambdaCoef(uint molA, uint box) const;

  // It's called in free energy calculation to calculate the change in
  // self energy in all lambda states
  virtual void ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const = 0;

  // It's called in free energy calculation to calculate the change in
  // correction energy in all lambda states
  virtual void ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                                const std::vector<double> &lambda_Coul,
                                const uint iState, const uint molIndex,
                                const uint box) const = 0;

  // It's called in free energy calculation to calculate the change in
  // reciprocal energy in all lambda states
  virtual void ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                           const std::vector<double> &lambda_Coul,
                           const uint iState, const uint molIndex,
                           const uint box) const = 0;



  bool electrostatic, ewald, multiParticleEnabled;
  double alpha;
  double recip_rcut, recip_rcut_Sq;

  std::vector<int> particleKind;
  std::vector<int> particleMol;

  // starting index of molecule
  std::vector<int> startMol;
  // length of each molecule
  std::vector<int> lengthMol;
  // starting atom index of each box
  int boxStart[BOX_TOTAL];
  // ending atom index of each box
  int boxEnd[BOX_TOTAL];
  // atom charges
  std::vector<double> particleCharge;
  // which atoms don't have charge
  std::vector<bool> particleHasNoCharge;
/* Other members */
};
#endif /*EWALD_H*/
