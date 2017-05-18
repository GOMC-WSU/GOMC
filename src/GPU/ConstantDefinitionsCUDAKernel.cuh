#ifndef CONSTANT_DEFINITIONS_CUDA_KERNEL
#define CONSTANT_DEFINITIONS_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "VariablesCUDA.cuh"

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
		       double const *epsilon_Cn, double const *n,
		       int VDW_Kind, int isMartini, int count,
		       double Rcut, double RcutLow, double Ron, double alpha,
		       int ewald, double diElectric_1);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
