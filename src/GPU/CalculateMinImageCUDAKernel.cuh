/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CALCULATE_MIN_IMAGE_CUDA_KERNEL_H
#define CALCULATE_MIN_IMAGE_CUDA_KERNEL_H

#ifdef GOMC_CUDA

#include "ConstantDefinitionsCUDAKernel.cuh"
#include <cuda.h>
#include <cuda_runtime.h>

__device__ inline double3 Difference3(const double *x, const double *y,
                                      const double *z, uint i, uint j) {
  return make_double3(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
}

__device__ inline void TransformSlantGPU(double3 &dist, const double3 &slant,
                                         const double *gpu_cell_x,
                                         const double *gpu_cell_y,
                                         const double *gpu_cell_z) {
  dist.x = slant.x * gpu_cell_x[0] + slant.y * gpu_cell_x[1] +
           slant.z * gpu_cell_x[2];
  dist.y = slant.x * gpu_cell_y[0] + slant.y * gpu_cell_y[1] +
           slant.z * gpu_cell_y[2];
  dist.z = slant.x * gpu_cell_z[0] + slant.y * gpu_cell_z[1] +
           slant.z * gpu_cell_z[2];
}

__device__ inline void TransformUnSlantGPU(double3 &dist, const double3 &slant,
                                           const double *gpu_Invcell_x,
                                           const double *gpu_Invcell_y,
                                           const double *gpu_Invcell_z) {
  dist.x = slant.x * gpu_Invcell_x[0] + slant.y * gpu_Invcell_x[1] +
           slant.z * gpu_Invcell_x[2];
  dist.y = slant.x * gpu_Invcell_y[0] + slant.y * gpu_Invcell_y[1] +
           slant.z * gpu_Invcell_y[2];
  dist.z = slant.x * gpu_Invcell_z[0] + slant.y * gpu_Invcell_z[1] +
           slant.z * gpu_Invcell_z[2];
}

__device__ inline void WrapPBC(double &v, const double &ax) {
  if (v >= ax)
    v -= ax;
  else if (v < 0)
    v += ax;
}

__device__ inline void WrapPBC3(double3 &v, const double3 &ax) {
  WrapPBC(v.x, ax.x);
  WrapPBC(v.y, ax.y);
  WrapPBC(v.z, ax.z);
}

__device__ inline void
WrapPBCNonOrth3(double3 &v, const double3 &ax, const double *gpu_cell_x,
                const double *gpu_cell_y, const double *gpu_cell_z,
                const double *gpu_Invcell_x, const double *gpu_Invcell_y,
                const double *gpu_Invcell_z) {
  double3 t;
  TransformUnSlantGPU(t, v, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  WrapPBC(t.x, ax.x);
  WrapPBC(t.y, ax.y);
  WrapPBC(t.z, ax.z);
  TransformSlantGPU(v, t, gpu_cell_x, gpu_cell_y, gpu_cell_z);
}

__device__ inline void UnwrapPBC(double &v, const double &ref, const double &ax,
                                 const double &halfax) {
  if (std::fabs(ref - v) > halfax) {
    if (ref < halfax)
      v -= ax;
    else
      v += ax;
  }
}

__device__ inline void UnwrapPBC3(double3 &v, const double3 &ref,
                                  const double3 &ax, const double3 &halfax) {
  UnwrapPBC(v.x, ref.x, ax.x, halfax.x);
  UnwrapPBC(v.y, ref.y, ax.y, halfax.y);
  UnwrapPBC(v.z, ref.z, ax.z, halfax.z);
}

__device__ inline void
UnwrapPBCNonOrth3(double3 &v, const double3 &ref, const double3 &ax,
                  const double3 &halfax, const double *gpu_cell_x,
                  const double *gpu_cell_y, const double *gpu_cell_z,
                  const double *gpu_Invcell_x, const double *gpu_Invcell_y,
                  const double *gpu_Invcell_z) {
  double3 t, tref;
  TransformUnSlantGPU(t, v, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  TransformUnSlantGPU(tref, ref, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  UnwrapPBC(t.x, tref.x, ax.x, halfax.x);
  UnwrapPBC(t.y, tref.y, ax.y, halfax.y);
  UnwrapPBC(t.z, tref.z, ax.z, halfax.z);
  TransformSlantGPU(v, t, gpu_cell_x, gpu_cell_y, gpu_cell_z);
}

__device__ inline double MinImageSignedGPU(double raw, const double ax,
                                           const double halfAx) {
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

__device__ inline double3 MinImageGPU(double3 rawVec, const double3 axis,
                                      const double3 halfAx) {
  rawVec.x = MinImageSignedGPU(rawVec.x, axis.x, halfAx.x);
  rawVec.y = MinImageSignedGPU(rawVec.y, axis.y, halfAx.y);
  rawVec.z = MinImageSignedGPU(rawVec.z, axis.z, halfAx.z);
  return rawVec;
}

__device__ inline double3
MinImageNonOrthGPU(double3 rawVec, const double3 &axis, const double3 &halfAx,
                   const double *gpu_cell_x, const double *gpu_cell_y,
                   const double *gpu_cell_z, const double *gpu_Invcell_x,
                   const double *gpu_Invcell_y, const double *gpu_Invcell_z) {
  double3 t;
  TransformUnSlantGPU(t, rawVec, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  t = MinImageGPU(t, axis, halfAx);
  TransformSlantGPU(rawVec, t, gpu_cell_x, gpu_cell_y, gpu_cell_z);
  return rawVec;
}

__device__ inline void
DeviceInRcut(double &distSq, double3 &dist, const double *gpu_x,
             const double *gpu_y, const double *gpu_z, int particleID,
             int otherParticle, double axx, double axy, double axz,
             int gpu_nonOrth, double *gpu_cell_x, double *gpu_cell_y,
             double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
             double *gpu_Invcell_z) {
  // calculate distance
  double3 axes, halfAx;
  dist.x = gpu_x[particleID] - gpu_x[otherParticle];
  dist.y = gpu_y[particleID] - gpu_y[otherParticle];
  dist.z = gpu_z[particleID] - gpu_z[otherParticle];

  axes.x = axx;
  halfAx.x = axx * 0.5;
  axes.y = axy;
  halfAx.y = axy * 0.5;
  axes.z = axz;
  halfAx.z = axz * 0.5;

  // minimum image
  if (gpu_nonOrth) {
    dist = MinImageNonOrthGPU(dist, axes, halfAx, gpu_cell_x, gpu_cell_y,
                              gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                              gpu_Invcell_z);
  } else {
    dist = MinImageGPU(dist, axes, halfAx);
  }

  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;
}

// Call by calculate energy whether it is in rCut
__device__ inline bool
InRcutGPU(double &distSq, const double *x, const double *y, const double *z,
          uint i, uint j, const double3 &axis, const double3 &halfAx,
          double gpu_rCut, int gpu_nonOrth, const double *gpu_cell_x,
          const double *gpu_cell_y, const double *gpu_cell_z,
          const double *gpu_Invcell_x, const double *gpu_Invcell_y,
          const double *gpu_Invcell_z) {
  double3 dist;
  dist = Difference3(x, y, z, i, j);
  // Do a binary print here of dist
  if (gpu_nonOrth) {
    dist = MinImageNonOrthGPU(dist, axis, halfAx, gpu_cell_x, gpu_cell_y,
                              gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                              gpu_Invcell_z);
  } else {
    dist = MinImageGPU(dist, axis, halfAx);
  }

  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;

  return ((gpu_rCut * gpu_rCut) > distSq);
}

// Call by force calculate to return the distance and virial component
__device__ inline bool
InRcutGPU(double &distSq, double3 &dist, const double *x, const double *y,
          const double *z, uint i, uint j, const double3 &axis,
          const double3 &halfAx, double gpu_rCut, int gpu_nonOrth,
          const double *gpu_cell_x, const double *gpu_cell_y,
          const double *gpu_cell_z, const double *gpu_Invcell_x,
          const double *gpu_Invcell_y, const double *gpu_Invcell_z) {
  dist = Difference3(x, y, z, i, j);
  if (gpu_nonOrth) {
    dist = MinImageNonOrthGPU(dist, axis, halfAx, gpu_cell_x, gpu_cell_y,
                              gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                              gpu_Invcell_z);
  } else {
    dist = MinImageGPU(dist, axis, halfAx);
  }

  distSq = dist.x * dist.x + dist.y * dist.y + dist.z * dist.z;

  return ((gpu_rCut * gpu_rCut) > distSq);
}

__device__ inline int FlatIndexGPU(int i, int j, int gpu_count) {
  return i + j * gpu_count;
}

__device__ inline double DotProductGPU(double kx, double ky, double kz,
                                       double x, double y, double z) {
  return (kx * x + ky * y + kz * z);
}

__device__ inline double DeviceGetLambdaVDW(int molA, int molB, int box,
                                            const bool *gpu_isFraction,
                                            const int *gpu_molIndex,
                                            const double *gpu_lambdaVDW) {
  double lambda = 1.0;
  if (gpu_isFraction[box]) {
    if (gpu_molIndex[box] == molA) {
      lambda *= gpu_lambdaVDW[box];
    }
    if (gpu_molIndex[box] == molB) {
      lambda *= gpu_lambdaVDW[box];
    }
  }
  return lambda;
}

__device__ inline double
DeviceGetLambdaCoulomb(int molA, int molB, int box, const bool *gpu_isFraction,
                       const int *gpu_molIndex,
                       const double *gpu_lambdaCoulomb) {
  double lambda = 1.0;
  if (gpu_isFraction[box]) {
    if (gpu_molIndex[box] == molA) {
      lambda *= gpu_lambdaCoulomb[box];
    }
    if (gpu_molIndex[box] == molB) {
      lambda *= gpu_lambdaCoulomb[box];
    }
  }
  return lambda;
}

__device__ inline double
DeviceGetLambdaCoulomb(int mol, int box, const bool *gpu_isFraction,
                       const int *gpu_molIndex,
                       const double *gpu_lambdaCoulomb) {
  double lambda = 1.0;
  if (gpu_isFraction[box]) {
    if (gpu_molIndex[box] == mol) {
      lambda = gpu_lambdaCoulomb[box];
    }
  }
  return lambda;
}

// Add atomic operations for GPUs that do not support it
// atomicAdd and atomicSub only support double for Compute Capability >= 6.0
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
static __inline__ __device__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
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
#endif /*CALCULATE_MIN_IMAGE_CUDA_KERNEL_H*/
