/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CONSTANT_DEFINITIONS_CUDA_KERNEL_H
#define CONSTANT_DEFINITIONS_CUDA_KERNEL_H

#ifdef GOMC_CUDA
#include "EnsemblePreprocessor.h"
#include "GeomLib.h"
#include "Molecules.h"
#include "VariablesCUDA.cuh"
#include <cuda.h>
#include <cuda_runtime.h>

const int GPU_VDW_STD_KIND = 0;
const int GPU_VDW_SHIFT_KIND = 1;
const int GPU_VDW_SWITCH_KIND = 2;
const int GPU_VDW_EXP6_KIND = 3;
const int MAX_PAIR_SIZE = 10000000;
const int NUMBER_OF_NEIGHBOR_CELLS = 27;
const int PARTICLES_PER_BLOCK = 32;
const int THREADS_PER_BLOCK = 128;
const int THREADS_PER_BLOCK_SM = 64;

void UpdateGPULambda(VariablesCUDA *vars, int *molIndex, double *lambdaVDW,
                     double *lambdaCoulomb, bool *isFraction);
void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n, int VDW_Kind,
                       int isMartini, int count, double Rcut, double RcutSq,
                       double const *rCutCoulomb, double const *rCutCoulombSq,
                       double RcutLow, double Ron, double const *alpha,
                       double const *alphaSq, int ewald, double diElectric_1);
void InitCoordinatesCUDA(VariablesCUDA *vars, uint maxAtomNumber,
                         uint maxAtomsInMol, uint maxMolNumber);
void InitEwaldVariablesCUDA(VariablesCUDA *vars, int numAtoms, uint imageTotal);
void InitExp6VariablesCUDA(VariablesCUDA *vars, double *rMin, double *expConst,
                           double *rMaxSq, uint size);
void InitMoleculeVariablesCUDA(VariablesCUDA *vars, const Molecules &mols);
void InitNeighborListVarsCUDA(VariablesCUDA *vars);
void InitPartVariablesCUDA(VariablesCUDA *vars,
                           const std::vector<int> &particleKind,
                           const std::vector<int> &particleMol,
                           const std::vector<double> &particleCharge);
void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void CopyRefToNewCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box);
void UpdateRecipCUDA(VariablesCUDA *vars, uint box);
void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z);
void UpdateInvCellBasisCUDA(VariablesCUDA *vars, uint box,
                            double *invCellBasis_x, double *invCellBasis_y,
                            double *invCellBasis_z);
void UpdateEnergyVecsCUDA(VariablesCUDA *vars, int newVecLen,
                          bool electrostatic);
void RebuildNeighborsCUDA(VariablesCUDA *vars,
                          std::vector<std::vector<int>> &neighbors,
                          const uint b);
void DestroyNeighborListCUDAVars(VariablesCUDA *vars);
void DestroyEwaldCUDAVars(VariablesCUDA *vars);
void DestroyExp6CUDAVars(VariablesCUDA *vars);
void DestroyCUDAVars(VariablesCUDA *vars);

#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL_H*/
