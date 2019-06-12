/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDAKernel.cuh"

__device__ inline void TransformSlantGPU(real &tx, real &ty, real &tz,
    real x, real y, real z,
    real *gpu_cell_x,
    real *gpu_cell_y,
    real *gpu_cell_z)
{
  tx = x * gpu_cell_x[0] + y * gpu_cell_x[1] + z * gpu_cell_x[2];
  ty = x * gpu_cell_y[0] + y * gpu_cell_y[1] + z * gpu_cell_y[2];
  tz = x * gpu_cell_z[0] + y * gpu_cell_z[1] + z * gpu_cell_z[2];
}

__device__ inline void TransformUnSlantGPU(real &tx, real &ty, real &tz,
    real x, real y, real z,
    real *gpu_Invcell_x,
    real *gpu_Invcell_y,
    real *gpu_Invcell_z)
{
  tx = x * gpu_Invcell_x[0] + y * gpu_Invcell_x[1] + z * gpu_Invcell_x[2];
  ty = x * gpu_Invcell_y[0] + y * gpu_Invcell_y[1] + z * gpu_Invcell_y[2];
  tz = x * gpu_Invcell_z[0] + y * gpu_Invcell_z[1] + z * gpu_Invcell_z[2];
}

__device__ inline real MinImageSignedGPU(real raw, real ax, real halfAx)
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

// Call by calculate energy whether it is in rCut
__device__ inline bool InRcutGPU(real &distSq, real gpu_x1, real gpu_y1,
                                 real gpu_z1, real gpu_x2, real gpu_y2,
                                 real gpu_z2, real xAxes, real yAxes,
                                 real zAxes, real xHalfAxes,
                                 real yHalfAxes, real zHalfAxes,
                                 real gpu_rCut, int gpu_nonOrth,
                                 real *gpu_cell_x, real *gpu_cell_y,
                                 real *gpu_cell_z, real *gpu_Invcell_x,
                                 real *gpu_Invcell_y, real *gpu_Invcell_z)
{
  distSq = 0;
  real tx, ty, tz;
  real dx = gpu_x1 - gpu_x2;
  real dy = gpu_y1 - gpu_y2;
  real dz = gpu_z1 - gpu_z2;

  if(gpu_nonOrth) {
    TransformUnSlantGPU(tx, ty, tz, dx, dy, dz, gpu_Invcell_x, gpu_Invcell_y,
                        gpu_Invcell_z);
    tx = MinImageSignedGPU(tx, xAxes, xHalfAxes);
    ty = MinImageSignedGPU(ty, yAxes, yHalfAxes);
    tz = MinImageSignedGPU(tz, zAxes, zHalfAxes);
    TransformSlantGPU(dx, dy, dz, tx, ty, tz, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    dx = MinImageSignedGPU(dx, xAxes, xHalfAxes);
    dy = MinImageSignedGPU(dy, yAxes, yHalfAxes);
    dz = MinImageSignedGPU(dz, zAxes, zHalfAxes);
  }

  distSq = dx * dx + dy * dy + dz * dz;
  return ((gpu_rCut * gpu_rCut) > distSq);
}

// Call by force calculate to return the distance and virial component
__device__ inline bool InRcutGPU(real &distSq, real &virX, real &virY,
                                 real &virZ, real gpu_x1, real gpu_y1,
                                 real gpu_z1, real gpu_x2, real gpu_y2,
                                 real gpu_z2, real xAxes, real yAxes,
                                 real zAxes, real xHalfAxes,
                                 real yHalfAxes, real zHalfAxes,
                                 real gpu_rCut, int gpu_nonOrth,
                                 real *gpu_cell_x, real *gpu_cell_y,
                                 real *gpu_cell_z, real *gpu_Invcell_x,
                                 real *gpu_Invcell_y, real *gpu_Invcell_z)
{
  distSq = 0;
  real tx, ty, tz;
  virX = gpu_x1 - gpu_x2;
  virY = gpu_y1 - gpu_y2;
  virZ = gpu_z1 - gpu_z2;
  if(gpu_nonOrth) {
    TransformUnSlantGPU(tx, ty, tz, virX, virY, virZ, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z);
    tx = MinImageSignedGPU(tx, xAxes, xHalfAxes);
    ty = MinImageSignedGPU(ty, yAxes, yHalfAxes);
    tz = MinImageSignedGPU(tz, zAxes, zHalfAxes);
    TransformSlantGPU(virX, virY, virZ, tx, ty, tz, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    virX = MinImageSignedGPU(virX, xAxes, xHalfAxes);
    virY = MinImageSignedGPU(virY, yAxes, yHalfAxes);
    virZ = MinImageSignedGPU(virZ, zAxes, zHalfAxes);
    distSq = virX * virX + virY * virY + virZ * virZ;
  }

  return ((gpu_rCut * gpu_rCut) > distSq);
}

__device__ inline int FlatIndexGPU(int i, int j, int gpu_count)
{
  return i + j * gpu_count;
}

__device__ inline real DotProductGPU(real kx, real ky, real kz,
                                       real x, real y, real z)
{
  return (kx * x + ky * y + kz * z);
}

// Add atomic operations for GPUs that do not support it
// atomicAdd and atomicSub only support double for Compute Capability >= 6.0
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
static __inline__ __device__ double atomicAdd(double *address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  if (val==0.0)
    return __longlong_as_double(old);
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val +__longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

#endif /*GOMC_CUDA*/
