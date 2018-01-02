/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
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
#define MAX_PAIR_SIZE 10000000

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n,
                       int VDW_Kind, int isMartini, int count,
                       double Rcut, double RcutLow, double Ron, double alpha,
                       int ewald, double diElectric_1);
void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
                         uint maxAtomsInMol, uint maxMolNumber);
void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal);
void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box);
void UpdateRecipCUDA(VariablesCUDA *vars, uint box);
void DestroyEwaldCUDAVars(VariablesCUDA *vars);
void DestroyCUDAVars(VariablesCUDA *vars);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
