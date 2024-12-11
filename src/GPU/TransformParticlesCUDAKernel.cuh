/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef TRANSFORM_PARTICLES_CUDA_KERNEL_H
#define TRANSFORM_PARTICLES_CUDA_KERNEL_H

#ifdef GOMC_CUDA
#include "Random123/philox.h"
#include <vector>
typedef r123::Philox4x64 RNG;

#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include "math.h"
#include <cuda.h>
#include <cuda_runtime.h>

void CallTranslateParticlesGPU(
    VariablesCUDA *vars, const std::vector<int8_t> &isMoleculeInvolved, int box,
    double t_max, double *mForcex, double *mForcey, double *mForcez,
    std::vector<int> &inForceRange, ulong step, unsigned int key, ulong seed,
    int atomCount, int molCount, double xAxes, double yAxes, double zAxes,
    XYZArray &newMolPos, XYZArray &newCOMs, double lambdaBETA, XYZArray &t_k,
    XYZArray &molForceRecRef);

void CallRotateParticlesGPU(
    VariablesCUDA *vars, const std::vector<int8_t> &isMoleculeInvolved, int box,
    double r_max, double *mTorquex, double *mTorquey, double *mTorquez,
    std::vector<int> &inForceRange, ulong step, unsigned int key, ulong seed,
    int atomCount, int molCount, double xAxes, double yAxes, double zAxes,
    XYZArray &newMolPos, XYZArray &newCOMs, double lambdaBETA, XYZArray &r_k);

__global__ void TranslateParticlesKernel(
    double t_max, double *molForcex, double *molForcey, double *molForcez,
    int *inForceRange, ulong step, unsigned int key, ulong seed, double *gpu_x,
    double *gpu_y, double *gpu_z, int *gpu_particleMol, int atomCount,
    double xAxes, double yAxes, double zAxes, double *gpu_comx,
    double *gpu_comy, double *gpu_comz, double *gpu_cell_x, double *gpu_cell_y,
    double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
    double *gpu_Invcell_z, int *gpu_nonOrth, double lambdaBETA,
    double *gpu_t_k_x, double *gpu_t_k_y, double *gpu_t_k_z,
    int8_t *gpu_isMoleculeInvolved, double *gpu_mForceRecx,
    double *gpu_mForceRecy, double *gpu_mForceRecz);

__global__ void RotateParticlesKernel(
    double r_max, double *molTorquex, double *molTorquey, double *molTorquez,
    int *inForceRange, ulong step, unsigned int key, ulong seed, double *gpu_x,
    double *gpu_y, double *gpu_z, int *gpu_particleMol, int atomCount,
    double xAxes, double yAxes, double zAxes, double *gpu_comx,
    double *gpu_comy, double *gpu_comz, double *gpu_cell_x, double *gpu_cell_y,
    double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
    double *gpu_Invcell_z, int *gpu_nonOrth, double lambdaBETA,
    double *gpu_r_k_x, double *gpu_r_k_y, double *gpu_r_k_z,
    int8_t *gpu_isMoleculeInvolved);

// Brownian Motion multiparticle
void BrownianMotionRotateParticlesGPU(
    VariablesCUDA *vars, const std::vector<unsigned int> &moleculeInvolved,
    XYZArray &mTorque, XYZArray &newMolPos, XYZArray &newCOMs, XYZArray &r_k,
    const XYZ &boxAxes, const double BETA, const double r_max, ulong step,
    unsigned int key, ulong seed, const int box, const bool isOrthogonal);

void BrownianMotionTranslateParticlesGPU(
    VariablesCUDA *vars, const std::vector<unsigned int> &moleculeInvolved,
    XYZArray &mForce, XYZArray &mForceRec, XYZArray &newMolPos,
    XYZArray &newCOMs, XYZArray &t_k, const XYZ &boxAxes, const double BETA,
    const double t_max, ulong step, unsigned int key, ulong seed, const int box,
    const bool isOrthogonal);

template <const bool isOrthogonal>
__global__ void BrownianMotionRotateKernel(
    int *startAtomIdx, double *gpu_x, double *gpu_y, double *gpu_z,
    double *molTorquex, double *molTorquey, double *molTorquez,
    double *gpu_comx, double *gpu_comy, double *gpu_comz, double *gpu_r_k_x,
    double *gpu_r_k_y, double *gpu_r_k_z, int *moleculeInvolved,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    double3 axis, double3 halfAx, int atomCount, double r_max, ulong step,
    unsigned int key, ulong seed, double BETA);

template <const bool isOrthogonal>
__global__ void BrownianMotionTranslateKernel(
    int *startAtomIdx, double *gpu_x, double *gpu_y, double *gpu_z,
    double *molForcex, double *molForcey, double *molForcez,
    double *molForceRecx, double *molForceRecy, double *molForceRecz,
    double *gpu_comx, double *gpu_comy, double *gpu_comz, double *gpu_t_k_x,
    double *gpu_t_k_y, double *gpu_t_k_z, int *moleculeInvolved,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    double3 axis, double3 halfAx, int atomCount, double t_max, ulong step,
    unsigned int key, ulong seed, double BETA);

#endif /*GOMC_CUDA*/
#endif /*TRANSFORM_PARTICLES_CUDA_KERNEL_H*/
