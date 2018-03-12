/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EWALDCACHED_H
#define EWALDCACHED_H

#include "BasicTypes.h"
#include "EnergyTypes.h"
#include <vector>
#include <stdio.h>
#include <cstring>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Molecules.h"
#include "Forcefield.h"
#include "Ewald.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "BoxDimensionsNonOrth.h"
#include "MoleculeKind.h"
#include "TrialMol.h"
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#endif

//
//    Calculating Electrostatic calculation with caching Fourier terms.
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

namespace cbmc
{
class TrialMol;
}
namespace config_setup
{
class SystemVals;
}

class EwaldCached : public Ewald
{
public:

  EwaldCached(StaticVals & stat, System & sys);
  ~EwaldCached();

  virtual void Init();

  virtual void AllocMem();

  //return size of image with defined Kmax value
  //uint GetImageSize();

  //initiliazie term used for ewald calculation
  void RecipInit(uint box, BoxDimensions const& boxAxes);

  //initiliazie wave vector for orthogonal box
  void RecipInitOrth(uint box, BoxDimensions const& boxAxes);

  //initiliazie wave vector for non-orthogonal box
  void RecipInitNonOrth(uint box, BoxDimensions const& boxAxes);

  //Get initial estimate of memory required
  void RecipCountInit(uint box, BoxDimensions const& boxAxes);

  //calculate self term for a box
  virtual double BoxSelf(BoxDimensions const& boxAxes, uint box) const;

  //setup reciprocate term for a box
  virtual void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

  //calculate reciprocate energy term for a box
  virtual double BoxReciprocal(uint box) const;

  //calculate reciprocate force term for a box
  virtual Virial ForceReciprocal(Virial& virial, uint box) const;


  //calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box)const;

  //calculate reciprocate term for displacement and rotation move
  virtual double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
                               const uint box, XYZ const*const newCOM = NULL);

  //calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol& trialMo) const;

  //calculate correction term after swap move
  virtual double SwapCorrection(const cbmc::TrialMol& trialMo) const;

  //calculate reciprocate term in destination box for swap move
  virtual double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int sourceBox, const int molIndex);

  //calculate reciprocate term in source box for swap move
  virtual double SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                 const uint box, const int molIndex);

  //back up reciptocate value to Ref (will be called during initialization)
  virtual void SetRecipRef(uint box);

  //update reciprocate values
  virtual void UpdateRecip(uint box);

  //update the hx,y,z hsqr and prefact
  virtual void UpdateRecipVec(uint box);

  //restore cosMol and sinMol
  virtual void RestoreMol(int molIndex);

  //update sinMol and cosMol
  virtual void exgMolCache();

  void UpdateVectorsAndRecipTerms();

private:
  
  double *cosMolRestore; //cos()*charge
  double *sinMolRestore; //sin()*charge
  double **cosMolRef;
  double **sinMolRef;
  double **cosMolBoxRecip;
  double **sinMolBoxRecip;
};


#endif /*EWALDCACHED_H*/
