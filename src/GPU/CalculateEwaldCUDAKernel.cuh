/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CALCULATE_EWALD_CUDA_KERNEL_H
#define CALCULATE_EWALD_CUDA_KERNEL_H

#ifdef GOMC_CUDA
#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars, XYZArray const &coords,
                               double const *kx, double const *ky,
                               double const *kz,
                               const std::vector<double> &molCharge,
                               uint imageSize, double *prefact, double *hsqr,
                               double &energyRecip, uint box);

void CallBoxReciprocalSumsGPU(VariablesCUDA *vars, XYZArray const &coords,
                              const std::vector<double> &molCharge,
                              uint imageSize, double &energyRecip, uint box);

void CallBoxForceReciprocalGPU(
    VariablesCUDA *vars, XYZArray &atomForceRec, XYZArray &molForceRec,
    const std::vector<double> &particleCharge,
    const std::vector<int> &particleMol, const std::vector<int> &particleUsed,
    double constValue, uint imageSize, XYZArray const &molCoords,
    BoxDimensions const &boxAxes, int moveType, int box);

void CallSwapReciprocalGPU(VariablesCUDA *vars, XYZArray const &coords,
                           const std::vector<double> &molCharge, uint imageSize,
                           double &energyRecip, uint box);

void CallMolReciprocalGPU(VariablesCUDA *vars, XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          const std::vector<double> &molCharge, uint imageSize,
                          double &energyRecip, uint box);

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars, uint imageSize, uint box,
                                  const std::vector<double> &molCharge,
                                  XYZArray const &molCoords,
                                  double &energyRecipNew, bool first_call);

// Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &molCharge,
                                      uint imageSize, double &energyRecip,
                                      const double lambdaCoef, uint box);

__global__ void BoxReciprocalSumsGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_molCharge, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_prefact,
    double *__restrict__ gpu_finalVal, int atomNumber);

__global__ void BoxForceReciprocalGPU(
    double *__restrict__ gpu_aForceRecx, double *__restrict__ gpu_aForceRecy,
    double *__restrict__ gpu_aForceRecz, double *__restrict__ gpu_mForceRecx,
    double *__restrict__ gpu_mForceRecy, double *__restrict__ gpu_mForceRecz,
    const double *__restrict__ gpu_particleCharge,
    const int *__restrict__ gpu_particleMol,
    const int *__restrict__ gpu_particleUsed,
    const int *__restrict__ gpu_startAtomIdx, const double *gpu_alpha,
    const double *__restrict__ gpu_alphaSq, const double constValue,
    const int imageSize, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_prefact,
    const double *__restrict__ gpu_sumRnew,
    const double *__restrict__ gpu_sumInew,
    const bool *__restrict__ gpu_isFraction,
    const int *__restrict__ gpu_molIndex,
    const double *__restrict__ gpu_lambdaCoulomb,
    const double *__restrict__ gpu_cell_x,
    const double *__restrict__ gpu_cell_y,
    const double *__restrict__ gpu_cell_z,
    const double *__restrict__ gpu_Invcell_x,
    const double *__restrict__ gpu_Invcell_y,
    const double *__restrict__ gpu_Invcell_z,
    const int *__restrict__ gpu_nonOrth, double axx, double axy, double axz,
    int moveType, int box);

__global__ void SwapReciprocalGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    int atomNumber, const double *__restrict__ gpu_molCharge,
    double *__restrict__ gpu_sumRnew, double *__restrict__ gpu_sumInew,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, int imageSize);

__global__ void MolReciprocalGPU(
    const double *__restrict__ gpu_cx, const double *__restrict__ gpu_cy,
    const double *__restrict__ gpu_cz, const double *__restrict__ gpu_nx,
    const double *__restrict__ gpu_ny, const double *__restrict__ gpu_nz,
    const double *__restrict__ gpu_kx, const double *__restrict__ gpu_ky,
    const double *__restrict__ gpu_kz, int atomNumber,
    const double *__restrict__ gpu_molCharge, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, int imageSize);

__global__ void MolExchangeReciprocalGPU(
    const double *__restrict__ gpu_kx, const double *__restrict__ gpu_ky,
    const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_prefactRef,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_molCharge,
    int numChargedParticles, const double *__restrict__ gpu_x,
    const double *__restrict__ gpu_y, const double *__restrict__ gpu_z,
    double *__restrict__ gpu_recipEnergies, int imageSize, bool first_call);

__global__ void ChangeLambdaMolReciprocalGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    int atomNumber, const double *__restrict__ gpu_molCharge,
    double *__restrict__ gpu_sumRnew, double *__restrict__ gpu_sumInew,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, double lambdaCoef, int imageSize);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL_H*/
