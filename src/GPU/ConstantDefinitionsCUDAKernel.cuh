/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
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
#define GPU_VDW_EXP6_KIND 3
#define MAX_PAIR_SIZE 10000000

void InitGPULambda(VariablesCUDA *vars, int *molIndex, int *kindIndex,
                   double *lambdaVDW, double *lambdaCoulomb, bool *isFraction);
void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n,
                       int VDW_Kind, int isMartini, int count,
                       double Rcut, double const *rCutCoulomb,
                       double RcutLow, double Ron, double const *alpha,
                       int ewald, double diElectric_1);
void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
                         uint maxAtomsInMol, uint maxMolNumber);
void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal);
void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box);
void UpdateRecipCUDA(VariablesCUDA *vars, uint box);
void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z);
void UpdateInvCellBasisCUDA(VariablesCUDA *vars, uint box,
                            double *invCellBasis_x, double *invCellBasis_y,
                            double *invCellBasis_z);
void DestroyEwaldCUDAVars(VariablesCUDA *vars);
void DestroyCUDAVars(VariablesCUDA *vars);
void InitExp6Variables(VariablesCUDA *vars, double *rMin, double *expConst,
                       double *rMaxSq, uint size);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
