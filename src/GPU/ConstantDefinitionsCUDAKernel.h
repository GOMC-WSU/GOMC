#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"

__constant__ double gpu_sigmaSq[1000];
__constant__ double gpu_epsilon_Cn[1000];
__constant__ double gpu_n[1000];
__constant__ int gpu_VDW_Kind;
__constant__ int gpu_count;
__constant__ bool gpu_isMartini;
__constant__ double gpu_rCut;
__constant__ double gpu_rCutLow;
__constant__ double gpu_rOn;
__constant__ double gpu_alpha;
__constant__ bool gpu_ewald;
__constant__ double gpu_diElectric_1;

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2

void InitGPUForceField(double const *sigmaSq, double const *epsilon_Cn,
		       double const *n, uint VDW_Kind, bool isMartini,
		       int count, int Rcut, int RcutLow, int Ron, double alpha,
		       bool ewald, double diElectric_1)
{
  int countSq = count * count;
  cudaMemcpyToSymbol("gpu_VDW_Kind", &VDW_Kind, sizeof(int));
  cudaMemcpyToSymbol("gpu_isMartini", &isMartini, sizeof(bool));
  cudaMemcpyToSymbol("gpu_sigmaSq", &sigmaSq, countSq * sizeof(double));
  cudaMemcpyToSymbol("gpu_epsilon_Cn", &epsilon_Cn, countSq * sizeof(double));
  cudaMemcpyToSymbol("gpu_n", &n, countSq * sizeof(double));
  cudaMemcpyToSymbol("gpu_rCut", &Rcut, sizeof(double));
  cudaMemcpyToSymbol("gpu_rCutLow", &RcutLow, sizeof(double));
  cudaMemcpyToSymbol("gpu_rOn", &Ron, sizeof(double));
  cudaMemcpyToSymbol("gpu_count", &count, sizeof(int));
  cudaMemcpyToSymbol("gpu_alpha", &alpha, sizeof(double));
  cudaMemcpyToSymbol("gpu_ewald", &ewald, sizeof(bool));
  cudaMemcpyToSymbol("gpu_diElectric_1", &diElectric_1, sizeof(double));
}

#endif
