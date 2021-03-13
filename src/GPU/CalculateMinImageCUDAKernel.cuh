/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDAKernel.cuh"

__device__ inline double3 Difference(double * x, double * y, double * z, uint i, uint j)
{
  return make_double3(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
}

__device__ inline void TransformSlantGPU(double3 & dist,
    double3 slant,
    double *gpu_cell_x,
    double *gpu_cell_y,
    double *gpu_cell_z)
{
  dist.x = slant.x * gpu_cell_x[0] + slant.y * gpu_cell_x[1] + slant.z * gpu_cell_x[2];
  dist.y = slant.x * gpu_cell_y[0] + slant.y * gpu_cell_y[1] + slant.z * gpu_cell_y[2];
  dist.z = slant.x * gpu_cell_z[0] + slant.y * gpu_cell_z[1] + slant.z * gpu_cell_z[2];
}

__device__ inline void TransformUnSlantGPU(double3 & dist,
    double3 slant,
    double *gpu_Invcell_x,
    double *gpu_Invcell_y,
    double *gpu_Invcell_z)
{
  dist.x = slant.x * gpu_Invcell_x[0] + slant.y * gpu_Invcell_x[1] + slant.z * gpu_Invcell_x[2];
  dist.y = slant.x * gpu_Invcell_y[0] + slant.y * gpu_Invcell_y[1] + slant.z * gpu_Invcell_y[2];
  dist.z = slant.x * gpu_Invcell_z[0] + slant.y * gpu_Invcell_z[1] + slant.z * gpu_Invcell_z[2];
}


__device__ inline void WrapPBC(double &v, double &ax)
{
  if(v >= ax)
    v -= ax;
  else if(v < 0)
    v += ax;
}

__device__ inline void WrapPBC_f3(double3 &v, double3 &ax)
{
  WrapPBC(v.x, ax.x);
  WrapPBC(v.y, ax.y);
  WrapPBC(v.z, ax.z);
}

__device__ inline void  UnwrapPBC(
  double &v,
  double &ref, 
  double &ax,
  double &halfax)
{
  if(abs(ref - v) > halfax) {
    if(ref < halfax)
      v -= ax;
    else
      v += ax;
  }
}

__device__ inline void UnwrapPBC_f3(
  double3 &v,
  double3 &ref, 
  double3 &ax,
  double3 &halfax)
{
  UnwrapPBC(v.x, ref.x, ax.x, halfax.x);
  UnwrapPBC(v.y, ref.y, ax.y, halfax.y);
  UnwrapPBC(v.z, ref.z, ax.z, halfax.z);
}

__device__ inline double MinImageSignedGPU(double raw, double ax, double halfAx)
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

__device__ inline void DeviceInRcut(
  double &distSq,
  double &distX,
  double &distY,
  double &distZ,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  int particleID,
  int otherParticle,
  double axx,
  double axy,
  double axz
)
{
  // calculate difference
  double diff_x = gpu_x[particleID] - gpu_x[otherParticle];
  double diff_y = gpu_y[particleID] - gpu_y[otherParticle];
  double diff_z = gpu_z[particleID] - gpu_z[otherParticle];

  // min image
  distX = MinImageSignedGPU(diff_x, axx, axx / 2.0);
  distY = MinImageSignedGPU(diff_y, axy, axy / 2.0);
  distZ = MinImageSignedGPU(diff_z, axz, axz / 2.0);

  distSq = distX * distX + distY * distY + distZ * distZ;
}

__device__ inline double3 MinImageGPU(double3 rawVec, double3 axis, double3 halfAx)
{
  rawVec.x = MinImageSignedGPU(rawVec.x, axis.x, halfAx.x);
  rawVec.y = MinImageSignedGPU(rawVec.y, axis.y, halfAx.y);
  rawVec.z = MinImageSignedGPU(rawVec.z, axis.z, halfAx.z);
  return rawVec;
}

// Call by calculate energy whether it is in rCut
__device__ inline bool InRcutGPU(double &distSq,
                                 double * x, double * y, double * z,
                                 uint i, uint j,
                                 double3 axis, double3 halfAx,
                                 double gpu_rCut, int gpu_nonOrth,
                                 double *gpu_cell_x, double *gpu_cell_y,
                                 double *gpu_cell_z, double *gpu_Invcell_x,
                                 double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double3 t, dist;
  t = make_double3(0.0, 0.0, 0.0);
  dist = make_double3(0.0, 0.0, 0.0);
  dist = Difference(x, y, z, i, j);
  // Do a binary print here of dist
  if(gpu_nonOrth) {
    TransformUnSlantGPU(t, dist, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z);
    t = MinImageGPU(t, axis, halfAx);
    TransformSlantGPU(dist, t, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    dist = MinImageGPU(dist, axis, halfAx);
  }

  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;

  return ((gpu_rCut * gpu_rCut) > distSq);
}

// Call by force calculate to return the distance and virial component
__device__ inline bool InRcutGPU(double &distSq, double3 & dist,
                                 double * x, double * y, double * z,
                                 uint i, uint j,
                                 double3 axis, double3 halfAx,
                                 double gpu_rCut, int gpu_nonOrth,
                                 double *gpu_cell_x, double *gpu_cell_y,
                                 double *gpu_cell_z, double *gpu_Invcell_x,
                                 double *gpu_Invcell_y, double *gpu_Invcell_z)
{
  distSq = 0;
  double3 t;
  t = make_double3(0.0, 0.0, 0.0);
  dist = Difference(x, y, z, i, j);
  // Do a binary print here of dist
  if(gpu_nonOrth) {
    TransformUnSlantGPU(t, dist, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z);
    t = MinImageGPU(t, axis, halfAx);
    TransformSlantGPU(dist, t, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z);
  } else {
    dist = MinImageGPU(dist, axis, halfAx);
  }

  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;

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

__device__ inline double DeviceGetLambdaCoulomb(int mol, int kind, int box,
    bool *gpu_isFraction,
    int *gpu_molIndex,
    int *gpu_kindIndex,
    double *gpu_lambdaCoulomb)
{
  double val = 1.0;
  if(gpu_isFraction[box]) {
    if((gpu_molIndex[box] == mol) && (gpu_kindIndex[box] == kind)) {
      val = gpu_lambdaCoulomb[box];
    }
  }
  return val;
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
