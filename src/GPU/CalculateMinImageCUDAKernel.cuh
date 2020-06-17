/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDAKernel.cuh"

__device__ inline void TransformSlantGPU(double &tx, double &ty, double &tz,
    double x, double y, double z,
    double *gpu_cell_x,
    double *gpu_cell_y,
    double *gpu_cell_z)
{
  tx = x * gpu_cell_x[0] + y * gpu_cell_x[1] + z * gpu_cell_x[2];
  ty = x * gpu_cell_y[0] + y * gpu_cell_y[1] + z * gpu_cell_y[2];
  tz = x * gpu_cell_z[0] + y * gpu_cell_z[1] + z * gpu_cell_z[2];
}

__device__ inline void TransformUnSlantGPU(double &tx, double &ty, double &tz,
    double x, double y, double z,
    double *gpu_Invcell_x,
    double *gpu_Invcell_y,
    double *gpu_Invcell_z)
{
  tx = x * gpu_Invcell_x[0] + y * gpu_Invcell_x[1] + z * gpu_Invcell_x[2];
  ty = x * gpu_Invcell_y[0] + y * gpu_Invcell_y[1] + z * gpu_Invcell_y[2];
  tz = x * gpu_Invcell_z[0] + y * gpu_Invcell_z[1] + z * gpu_Invcell_z[2];
}

__device__ inline double MinImageSignedGPU(double raw, double ax, double halfAx)
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

// Call by calculate energy whether it is in rCut
__device__ inline bool InRcutGPU(double &distSq, double gpu_x1, double gpu_y1,
                                 double gpu_z1, double gpu_x2, double gpu_y2,
                                 double gpu_z2, double xAxes, double yAxes,
                                 double zAxes, double xHalfAxes,
                                 double yHalfAxes, double zHalfAxes,
                                 double gpu_rCut, int gpu_nonOrth,
                                 double *gpu_cell_x, double *gpu_cell_y,
                                 double *gpu_cell_z, double *gpu_Invcell_x,
                                 double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double tx, ty, tz;
  double dx = gpu_x1 - gpu_x2;
  double dy = gpu_y1 - gpu_y2;
  double dz = gpu_z1 - gpu_z2;

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
__device__ inline bool InRcutGPU(double &distSq, double &virX, double &virY,
                                 double &virZ, double gpu_x1, double gpu_y1,
                                 double gpu_z1, double gpu_x2, double gpu_y2,
                                 double gpu_z2, double xAxes, double yAxes,
                                 double zAxes, double xHalfAxes,
                                 double yHalfAxes, double zHalfAxes,
                                 double gpu_rCut, int gpu_nonOrth,
                                 double *gpu_cell_x, double *gpu_cell_y,
                                 double *gpu_cell_z, double *gpu_Invcell_x,
                                 double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double tx, ty, tz;
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

__device__ inline double DotProductGPU(double kx, double ky, double kz,
                                       double x, double y, double z)
{
  return (kx * x + ky * y + kz * z);
}

__device__ inline double DeviceGetLambdaVDW(int molA, int kindA, int molB,
    int kindB, int box,
    bool *gpu_isFraction,
    int *gpu_molIndex,
    int *gpu_kindIndex,
    double *gpu_lambdaVDW)
{
  double lambda = 1.0;
  if(gpu_isFraction[box]) {
    if((gpu_molIndex[box] == molA) && (gpu_kindIndex[box] == kindA)) {
      lambda *= gpu_lambdaVDW[box];
    }
    if((gpu_molIndex[box] == molB) && (gpu_kindIndex[box] == kindB)) {
      lambda *= gpu_lambdaVDW[box];
    }
  }
  return lambda;
}

__device__ inline double DeviceGetLambdaCoulomb(int molA, int kindA, int molB,
    int kindB, int box,
    bool *gpu_isFraction,
    int *gpu_molIndex,
    int *gpu_kindIndex,
    double *gpu_lambdaCoulomb)
{
  double lambda = 1.0;
  if(gpu_isFraction[box]) {
    if((gpu_molIndex[box] == molA) && (gpu_kindIndex[box] == kindA)) {
      lambda *= gpu_lambdaCoulomb[box];
    }
    if((gpu_molIndex[box] == molB) && (gpu_kindIndex[box] == kindB)) {
      lambda *= gpu_lambdaCoulomb[box];
    }
  }
  return lambda;
}

// Add atomic operations for GPUs that do not support it
// atomicAdd and atomicSub only support double for Compute Capability >= 6.0
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
static __inline__ __device__ double atomicAdd(double *address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  if (val == 0.0)
    return __longlong_as_double(old);
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

#endif /*GOMC_CUDA*/
