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

  void Init();

  void AllocMem();

  //setup reciprocate term for a box
  void BoxReciprocalSetup(uint box, XYZArray const& molCoords);

  //calculate reciprocate energy term for a box
  double BoxReciprocal(uint box) const;

  //calculate reciprocate term for displacement and rotation move
  double MolReciprocal(XYZArray const& molCoords, const uint molIndex,
                               const uint box, XYZ const*const newCOM = NULL);

  //calculate reciprocate term in destination box for swap move
  double SwapDestRecip(const cbmc::TrialMol &newMol, const uint box,
                               const int sourceBox, const int molIndex);

  //calculate reciprocate term in source box for swap move
  double SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                 const uint box, const int molIndex);

  //restore cosMol and sinMol
  void RestoreMol(int molIndex);

  //update sinMol and cosMol
  void exgMolCache();

private:
  
  double *cosMolRestore; //cos()*charge
  double *sinMolRestore; //sin()*charge
  double **cosMolRef;
  double **sinMolRef;
  double **cosMolBoxRecip;
  double **sinMolBoxRecip;
};


#endif /*EWALDCACHED_H*/
