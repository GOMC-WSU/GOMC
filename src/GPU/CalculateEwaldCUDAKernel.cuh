/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
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

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars,
                               XYZArray const &coords,
                               double const *kx,
                               double const *ky,
                               double const *kz,
                               std::vector<double> particleCharge,
                               uint imageSize,
                               double *sumRnew,
                               double *sumInew,
                               double *prefact,
                               double *hsqr,
                               double &energyRecip,
                               uint box);

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          std::vector<double> particleCharge,
                          uint imageSize,
                          double *sumRnew,
                          double *sumInew,
                          double &energyRecipNew,
                          uint box);

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           std::vector<double> particleCharge,
                           uint imageSize,
                           double *sumRnew,
                           double *sumInew,
                           int const insert,
                           double &energyRecipNew,
                           uint box);

__global__ void BoxReciprocalSetupGPU(double * gpu_x,
                                      double * gpu_y,
                                      double * gpu_z,
                                      double * gpu_kx,
                                      double * gpu_ky,
                                      double * gpu_kz,
                                      double atomNumber,
                                      double * gpu_particleCharge,
                                      double * gpu_sumRnew,
                                      double * gpu_sumInew,
                                      int imageSize);

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber,
                                 double *gpu_particleCharge,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_sumRref,
                                 double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_energyRecipNew,
                                 int imageSize);

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                  double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                  int atomNumber,
                                  double *gpu_particleCharge,
                                  double *gpu_sumRnew,
                                  double *gpu_sumInew,
                                  double *gpu_sumRref,
                                  double *gpu_sumIref,
                                  double *gpu_prefactRef,
                                  int insert,
                                  double *gpu_energyRecipNew,
                                  int imageSize);

__global__ void BoxReciprocalGPU(double *gpu_prefact,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_energyRecip,
                                 int imageSize);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL*/
