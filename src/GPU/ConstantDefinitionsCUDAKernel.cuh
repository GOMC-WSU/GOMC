/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CONSTANT_DEFINITIONS_CUDA_KERNEL
#define CONSTANT_DEFINITIONS_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "VariablesCUDA.cuh"
#include "EnsemblePreprocessor.h"

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2
#define GPU_VDW_EXP6_KIND 3
#define MAX_PAIR_SIZE 10000000

void UpdateGPULambda(VariablesCUDA *vars, int *molIndex, double *lambdaVDW,
                    double *lambdaCoulomb, bool *isFraction);
void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n,
                       int VDW_Kind, int isMartini, int count,
                       double Rcut, double const *rCutCoulomb,
                       double RcutLow, double Ron, double const *alpha,
                       int ewald, double diElectric_1);
void InitCoordinatesCUDA(VariablesCUDA *vars, uint maxAtomNumber,
                         uint maxAtomsInMol, uint maxMolNumber);
void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal);
void InitExp6VariablesCUDA(VariablesCUDA *vars, double *rMin, double *expConst,
                           double *rMaxSq, uint size);
void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void CopyRefToNewCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box);
void UpdateRecipCUDA(VariablesCUDA *vars, uint box);
void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z);
void UpdateInvCellBasisCUDA(VariablesCUDA *vars, uint box,
                            double *invCellBasis_x, double *invCellBasis_y,
                            double *invCellBasis_z);
void DestroyEwaldCUDAVars(VariablesCUDA *vars);
void DestroyExp6CUDAVars(VariablesCUDA *vars);
void DestroyCUDAVars(VariablesCUDA *vars);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
