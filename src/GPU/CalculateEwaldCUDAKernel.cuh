/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CALCULATE_EWALD_CUDA_KERNEL
#define CALCULATE_EWALD_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include <vector>
#include "XYZArray.h"

using namespace std;

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars,
                               XYZArray const &coords,
                               real const *kx,
                               real const *ky,
                               real const *kz,
                               vector<real> particleCharge,
                               uint imageSize,
                               real *sumRnew,
                               real *sumInew,
                               real *prefact,
                               real *hsqr,
                               real &energyRecip,
                               uint box);

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          vector<real> particleCharge,
                          uint imageSize,
                          real *sumRnew,
                          real *sumInew,
                          real &energyRecipNew,
                          uint box);

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           vector<real> particleCharge,
                           uint imageSize,
                           real *sumRnew,
                           real *sumInew,
                           int const insert,
                           real &energyRecipNew,
                           uint box);

__global__ void BoxReciprocalSetupGPU(real * gpu_x,
                                      real * gpu_y,
                                      real * gpu_z,
                                      real * gpu_kx,
                                      real * gpu_ky,
                                      real * gpu_kz,
                                      real atomNumber,
                                      real * gpu_particleCharge,
                                      real * gpu_sumRnew,
                                      real * gpu_sumInew,
                                      int imageSize);

__global__ void MolReciprocalGPU(real *gpu_cx, real *gpu_cy, real *gpu_cz,
                                 real *gpu_nx, real *gpu_ny, real *gpu_nz,
                                 real *gpu_kx, real *gpu_ky, real *gpu_kz,
                                 int atomNumber,
                                 real *gpu_particleCharge,
                                 real *gpu_sumRnew,
                                 real *gpu_sumInew,
                                 real *gpu_sumRref,
                                 real *gpu_sumIref,
                                 real *gpu_prefactRef,
                                 real *gpu_energyRecipNew,
                                 int imageSize);

__global__ void SwapReciprocalGPU(real *gpu_x, real *gpu_y, real *gpu_z,
                                  real *gpu_kx, real *gpu_ky, real *gpu_kz,
                                  int atomNumber,
                                  real *gpu_particleCharge,
                                  real *gpu_sumRnew,
                                  real *gpu_sumInew,
                                  real *gpu_sumRref,
                                  real *gpu_sumIref,
                                  real *gpu_prefactRef,
                                  int insert,
                                  real *gpu_energyRecipNew,
                                  int imageSize);

__global__ void BoxReciprocalGPU(real *gpu_prefact,
                                 real *gpu_sumRnew,
                                 real *gpu_sumInew,
                                 real *gpu_energyRecip,
                                 int imageSize);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL*/
