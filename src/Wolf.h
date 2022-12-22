/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef WOLF_H
#define WOLF_H

#include <stdio.h>

#include <cassert>
#include <cstring>
#include <vector>

#include "ElectrostaticBase.h"


#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "Forcefield.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "TrialMol.h"
#ifdef _OPENMP
#include <omp.h>
#endif
//
//    Calculating Wolf Approximation of the Electrostatic
//    functions for self and correction terms.
//
//    Developed by G. Schwing and Mohammad S. Barhaghi
//
//
typedef double BOX_SIZE_DOUBLE_ARRAY[BOX_TOTAL];

class StaticVals;
class System;
class Forcefield;
class Molecules;
class MoleculeLookup;
class MoleculeKind;
class Coordinates;
class COM;
class XYZArray;
class BoxDimensions;
class CalculateEnergy;
class Lambda;

enum wolfKind { WOLF_RAHBARI_KIND, WOLF_WAIBEL2018_KIND, WOLF_WAIBEL2019_KIND };

class Wolf : public ElectrostaticBase{
  // friend class CalculateEnergy;
public:
  Wolf(StaticVals &stat, System &sys);
  virtual ~Wolf() {}

  virtual void Init();

  virtual void AllocMem() {}

  // initialize term used for ewald calculation
  virtual void RecipInit(uint box, BoxDimensions const &boxAxes) {}

  // initialize wave vector for orthogonal box
  virtual void RecipInitOrth(uint box, BoxDimensions const &boxAxes) {}

  // initialize wave vector for non-orthogonal box
  virtual void RecipInitNonOrth(uint box, BoxDimensions const &boxAxes) {}

  // Get initial estimate of memory required
  virtual void RecipCountInit(uint box, BoxDimensions const &boxAxes) {}

  // compute reciprocal term for a box with a new volume
  virtual void BoxReciprocalSetup(uint box, XYZArray const &molCoords) {}

  // compute reciprocal term for a box when not testing a volume change
  virtual void BoxReciprocalSums(uint box, XYZArray const &molCoords) {}

  // calculate reciprocal energy term for a box
  virtual double BoxReciprocal(uint box, bool isNewVolume) const { return 0.0; }

  // calculate self term for a box
  virtual double BoxSelf(uint box) const;

  // calculate reciprocal force term for a box
  virtual Virial VirialReciprocal(Virial &virial, uint box) const { return Virial(); }

  // calculate reciprocal term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const &molCoords, const uint molIndex,
                               const uint box) { return 0.0; }

  // calculate reciprocal term for lambdaNew and Old with same coordinates
  virtual double ChangeLambdaRecip(XYZArray const &molCoords,
                                   const double lambdaOld,
                                   const double lambdaNew, const uint molIndex,
                                   const uint box) { return 0.0; }

  // calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box) const;

  // calculate reciprocal term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex) { return 0.0; }

  // calculate reciprocal term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol, const uint box,
                                 const int molIndex) { return 0.0; }

  // calculate reciprocal term for inserting some molecules (kindA) in
  // destination box and removing a molecule (kindB) from destination box
  virtual double
  MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                        const std::vector<cbmc::TrialMol> &oldMol,
                        const std::vector<uint> &molIndexNew,
                        const std::vector<uint> &molIndexOld, bool first_call) { return 0.0; }

  // calculate correction term after swap move, with lambda = 1
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol) const;

  // calculate correction term after swap move, with system lambda
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol,
                                const uint molIndex) const;

  // back up reciprocal value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box) {}

  // update reciprocal values
  virtual void UpdateRecip(uint box) {}

  // copy reciprocal values from ref to new
  virtual void CopyRecip(uint box) {}

  // update kx, ky, kz, hsqr and prefact
  virtual void UpdateRecipVec(uint box) {}

  // calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol &trialMol) const;

  // restore cosMol and sinMol
  virtual void RestoreMol(int molIndex) {}

  // Find the largest kvector
  virtual uint findLargeImage() { return 0;}

  // update sinMol and cosMol
  virtual void exgMolCache() {}

  // backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache() {}

  /// This function performs three actions:
  /// 1. Initialize k vectors
  /// 2. Run BoxReciprocalSetup to calculate sumRnew and sumInew vectors
  /// 3. Then copy all new vectors to ref vectors
  /// @param bool: whether to print the vector size
  virtual void UpdateVectorsAndRecipTerms(bool output) {}

  // calculate reciprocal force term for a box with molCoords
  virtual void BoxForceReciprocal(XYZArray const &molCoords,
                                  XYZArray &atomForceRec, XYZArray &molForceRec,
                                  uint box) {}

  double GetLambdaCoef(uint molA, uint box) const;

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
                           const uint box) const {}

private:
  const uint &wolfKind, &coulKind;
  const BOX_SIZE_DOUBLE_ARRAY& wolfAlpha; //alpha term for Wolf Electrostatic and constant factors
  const BOX_SIZE_DOUBLE_ARRAY& wolfFactor1; //alpha term for Wolf Electrostatic and constant factors
  const BOX_SIZE_DOUBLE_ARRAY& wolfFactor2;  //alpha term for Wolf Electrostatic and constant factors
  const BOX_SIZE_DOUBLE_ARRAY& wolfFactor3; //alpha term for Wolf Electrostatic and constant factors  
  const BOX_SIZE_DOUBLE_ARRAY& rCutCoulomb;  //alpha term for Wolf Electrostatic and constant factors
  const BOX_SIZE_DOUBLE_ARRAY& rCutCoulombSq; //alpha term for Wolf Electrostatic and constant factors  
  bool oneThree, oneFour;
  double scaling_14;
  double RahbariCorrection(MoleculeKind &thisKind,
                                uint box,
                                XYZ virComponents,
                                uint atomSize,
                                uint start) const;
  double Waibel2018Correction(MoleculeKind &thisKind,
                                uint box,
                                XYZ virComponents,
                                uint atomSize,
                                uint start) const;
  double Waibel2019Correction(MoleculeKind &thisKind,
                                uint box,
                                XYZ virComponents,
                                uint atomSize,
                                uint start) const;

  double RahbariCorrection(const cbmc::TrialMol &trialMol) const;
  double Waibel2018Correction(const cbmc::TrialMol &trialMol) const;
  double Waibel2019Correction(const cbmc::TrialMol &trialMol) const;

protected:
  const Forcefield &ff;
  const Molecules &mols;
  const Coordinates &currentCoords;
  const MoleculeLookup &molLookup;
  const BoxDimensions &currentAxes;
  const COM &currentCOM;
  const SystemPotential &sysPotRef;
  const Lambda &lambdaRef;

};

#endif /*EWALD_H*/
