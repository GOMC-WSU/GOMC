/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
//    Calculating self, correction and reciprocate part of ewald
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


class Ewald
{
  //friend class CalculateEnergy;
public:

  Ewald(StaticVals & stat, System & sys);
  ~Ewald();

  virtual void Init();

  virtual void AllocMem();

  //initiliazie term used for ewald calculation
  virtual void RecipInit(uint box, BoxDimensions const& boxAxes);

  //initiliazie wave vector for orthogonal box
  virtual void RecipInitOrth(uint box, BoxDimensions const& boxAxes);

  //initiliazie wave vector for non-orthogonal box
  virtual void RecipInitNonOrth(uint box, BoxDimensions const& boxAxes);

  //Get initial estimate of memory required
  void RecipCountInit(uint box, BoxDimensions const& boxAxes);

  //setup reciprocate term for a box
  virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

  //calculate reciprocate energy term for a box
  virtual double BoxReciprocal(uint box) const;

  //calculate self term for a box
  virtual double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

  //calculate reciprocate force term for a box
  virtual Virial ForceReciprocal(Virial& virial, uint box) const;

  //calculate reciprocate term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
                               const uint box);

  //calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box)const;

  //calculate reciprocate term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int molIndex);

  //calculate reciprocate term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                 const uint box, const int molIndex);

  //calculate reciprocate term for inserting some molecules (kindA) in
  //destination box and removing a molecule (kindB) from destination box
  virtual double SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
                           const std::vector<cbmc::TrialMol> &oldMol);

  //calculate correction term after swap move
  virtual double SwapCorrection(const cbmc::TrialMol& trialMol) const;

  //back up reciptocate value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box);

  //update reciprocate values
  virtual void UpdateRecip(uint box);

  //update the hx,y,z hsqr and prefact
  virtual void UpdateRecipVec(uint box);

  //calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol& trialMo) const;

  //restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  //Find the largest kvector
  uint findLargeImage();

  //update sinMol and cosMol
  virtual void exgMolCache();

  //backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
  virtual void backupMolCache();

  virtual void UpdateVectorsAndRecipTerms();

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

  uint *imageSize;
  uint *imageSizeRef;
  //const uint imageTotal = GetImageSize();
  uint imageTotal;
  uint *kmax;
  double **sumRnew; //cosine serries
  double **sumInew; //sine serries
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
};



#endif /*EWALD_H*/
