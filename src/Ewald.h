/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EWALD_H
#define EWALD_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include "Molecules.h"
#include "Forcefield.h"
#include "TrialMol.h"
#include "MoleculeLookup.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif
//
//    Calculating Electrostatic calculation without caching Fourier terms.
//    Energy Calculation functions for Ewald summation method
//    Calculating self, correction and reciprocal part of ewald
//
//    Developed by Y. Li and Mohammad S. Barhaghi
//
//

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


class Ewald
{
  //friend class CalculateEnergy;
public:

  Ewald(StaticVals & stat, System & sys);
  virtual ~Ewald();

  virtual void Init();

  virtual void AllocMem();

  //initialize term used for ewald calculation
  virtual void RecipInit(uint box, BoxDimensions const& boxAxes);

  //initialize wave vector for orthogonal box
  virtual void RecipInitOrth(uint box, BoxDimensions const& boxAxes);

  //initialize wave vector for non-orthogonal box
  virtual void RecipInitNonOrth(uint box, BoxDimensions const& boxAxes);

  //Get initial estimate of memory required
  void RecipCountInit(uint box, BoxDimensions const& boxAxes);

  //setup reciprocal term for a box
  virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

  //calculate reciprocal energy term for a box
  virtual double BoxReciprocal(uint box) const;

  //calculate self term for a box
  virtual double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

  //calculate reciprocal force term for a box
  virtual Virial VirialReciprocal(Virial& virial, uint box) const;

  //calculate reciprocal term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
                               const uint box);

  //calculate reciprocal term for lambdaNew and Old with same coordinates
  virtual double CFCMCRecip(XYZArray const& molCoords, const double lambdaOld,
                            const double lambdaNew, const uint molIndex,
                            const uint box);

  //calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box)const;

  //calculate reciprocal term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex);

  //calculate reciprocal term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                 const uint box, const int molIndex);

  //calculate reciprocal term for inserting some molecules (kindA) in
  //destination box and removing a molecule (kindB) from destination box
  virtual double SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
                           const std::vector<cbmc::TrialMol> &oldMol,
                           const std::vector<uint> molIndexNew,
                           const std::vector<uint> molIndexOld);

  //calculate correction term after swap move, with lambda = 1
  virtual double SwapCorrection(const cbmc::TrialMol& trialMol) const;

  //calculate correction term after swap move, with system lambda
  virtual double SwapCorrection(const cbmc::TrialMol& trialMol,
                                const uint molIndex) const;

  //back up reciprocal value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box);

  //update reciprocal values
  virtual void UpdateRecip(uint box);

  //update the hx,y,z hsqr and prefact
  virtual void UpdateRecipVec(uint box);

  //calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol& trialMol) const;

  //restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  //Find the largest kvector
  uint findLargeImage();

  //update sinMol and cosMol
  virtual void exgMolCache();

  //backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache();

  virtual void UpdateVectorsAndRecipTerms();

  //calculate reciprocal force term for a box with molCoords
  virtual void BoxForceReciprocal(XYZArray const& molCoords,
                                  XYZArray& atomForceRec,
                                  XYZArray& molForceRec,
                                  uint box);

  double GetLambdaCoef(uint molA, uint box) const;

  //It's called in free energy calculation to calculate the change in
  // self energy in all lambda states
  virtual void ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const;

  //It's called in free energy calculation to calculate the change in
  // correction energy in all lambda states
  virtual void ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                                const std::vector<double> &lambda_Coul,
                                const uint iState, const uint molIndex,
                                const uint box) const;

  //It's called in free energy calculation to calculate the change in
  // reciprocal energy in all lambda states
  virtual void ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                           const std::vector<double> &lambda_Coul,
                           const uint iState, const uint molIndex,
                           const uint box) const;

private:
  double currentEnergyRecip[BOXES_WITH_U_NB];

protected:
  const Forcefield& ff;
  const Molecules& mols;
  const Coordinates& currentCoords;
  const MoleculeLookup& molLookup;
  const BoxDimensions& currentAxes;
  const COM& currentCOM;
  const SystemPotential &sysPotRef;
  const Lambda& lambdaRef;

  bool electrostatic, ewald, multiParticleEnabled;
  double alpha;
  double recip_rcut, recip_rcut_Sq;
  uint *imageSize;
  uint *imageSizeRef;
  //const uint imageTotal = GetImageSize();
  uint imageTotal;
  uint *kmax;
  double **sumRnew; //cosine series
  double **sumInew; //sine series
  double **sumRref;
  double **sumIref;

  double **kx, **kxRef;
  double **ky, **kyRef;
  double **kz, **kzRef;
  double **hsqr, **hsqrRef;
  double **prefact, **prefactRef;


  std::vector<int> particleKind;
  std::vector<int> particleMol;
  std::vector<double> particleCharge;
  std::vector<bool> particleHasNoCharge;

};



#endif /*EWALD_H*/
