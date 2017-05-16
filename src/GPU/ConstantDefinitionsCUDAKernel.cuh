#ifndef CONSTANT_DEFINITIONS_CUDA_KERNEL
#define CONSTANT_DEFINITIONS_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "VariablesCUDA.h"

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
		       double const *epsilon_Cn, double const *n,
		       uint VDW_Kind, bool isMartini, int count,
		       int Rcut, int RcutLow, int Ron, double alpha,
		       bool ewald, double diElectric_1);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
