#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include <iostream>
#include "VariablesCUDA.h"

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2

void InitGPUForceField(VariablesCUDA& vars, double const *sigmaSq, 
		       double const *epsilon_Cn,
		       double const *n, int VDW_Kind, int isMartini,
		       int count, double Rcut, double RcutLow,
		       double Ron, double alpha,
		       int ewald, double diElectric_1)
{
  int countSq = count * count;
  cudaMalloc(&vars.gpu_sigmaSq, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_epsilon_Cn, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_n, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_VDW_Kind, sizeof(int));
  cudaMalloc(&vars.gpu_isMartini, sizeof(int));
  cudaMalloc(&vars.gpu_count, sizeof(int));
  cudaMalloc(&vars.gpu_rCut, sizeof(double));
  cudaMalloc(&vars.gpu_rCutLow, sizeof(double));
  cudaMalloc(&vars.gpu_rOn, sizeof(double));
  cudaMalloc(&vars.gpu_alpha, sizeof(double));
  cudaMalloc(&vars.gpu_ewald, sizeof(int));
  cudaMalloc(&vars.gpu_diElectric_1, sizeof(double));

  cudaMemcpy(vars.gpu_sigmaSq, sigmaSq, countSq * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_epsilon_Cn, epsilon_Cn, countSq * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_n, n, countSq * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_VDW_Kind, &VDW_Kind, sizeof(int),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_isMartini, &isMartini, sizeof(int),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_count, &count, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCut, &Rcut, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutLow, &RcutLow, sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rOn, &Ron, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_alpha, &alpha, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_ewald, &ewald, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_diElectric_1, &diElectric_1, sizeof(double),
	     cudaMemcpyHostToDevice);
}

#endif /*GOMC_CUDA*/
