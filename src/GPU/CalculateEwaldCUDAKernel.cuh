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

void CallBoxForceReciprocalGPU(
    VariablesCUDA *vars, XYZArray &atomForceRec, XYZArray &molForceRec,
    const std::vector<double> &particleCharge,
    const std::vector<int> &particleMol, const std::vector<int> &particleUsed,
    double constValue, uint imageSize, XYZArray const &molCoords,
    BoxDimensions const &boxAxes, int moveType, int box);

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars, XYZArray const &coords,
                               double const *kx, double const *ky,
                               double const *kz,
                               const std::vector<double> &molCharge,
                               uint imageSize, double *prefact, double *hsqr,
                               double &energyRecip, uint box);

void CallBoxReciprocalSumsGPU(VariablesCUDA *vars, XYZArray const &coords,
                              const std::vector<double> &molCharge,
                              uint imageSize, double &energyRecip, uint box);

void CallMolReciprocalGPU(VariablesCUDA *vars, XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          const std::vector<double> &molCharge, uint imageSize,
                          double &energyRecip, uint box);

// Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &molCharge,
                                      uint imageSize, double &energyRecip,
                                      const double lambdaCoef, uint box);

void CallSwapReciprocalGPU(VariablesCUDA *vars, XYZArray const &coords,
                           const std::vector<double> &molCharge, uint imageSize,
                           double &energyRecip, uint box);

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars, uint imageSize, uint box,
                                  const std::vector<double> &molCharge,
                                  XYZArray const &molCoords,
                                  double &energyRecipNew, bool first_call);

__global__ void BoxForceReciprocalGPU(
    double *gpu_aForceRecx, double *gpu_aForceRecy, double *gpu_aForceRecz,
    double *gpu_mForceRecx, double *gpu_mForceRecy, double *gpu_mForceRecz,
    const double *gpu_particleCharge, const int *gpu_particleMol,
    const int *gpu_particleUsed, const int *gpu_startAtomIdx, double *gpu_alpha,
    double *gpu_alphaSq, double constValue, int imageSize, double *gpu_kx,
    double *gpu_ky, double *gpu_kz, double *gpu_x, double *gpu_y, double *gpu_z,
    double *gpu_prefact, double *gpu_sumRnew, double *gpu_sumInew,
    bool *gpu_isFraction, int *gpu_molIndex, double *gpu_lambdaCoulomb,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    int *gpu_nonOrth, double axx, double axy, double axz, int moveType,
    int box);

__global__ void BoxReciprocalSumsGPU(double *gpu_x, double *gpu_y,
                                     double *gpu_z, double *gpu_kx,
                                     double *gpu_ky, double *gpu_kz,
                                     double *gpu_molCharge, double *gpu_sumRnew,
                                     double *gpu_sumInew,
                                     double *gpu_prefactRef,
                                     double *gpu_finalVal, int atomNumber);

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber, double *gpu_molCharge,
                                 double *gpu_sumRnew, double *gpu_sumInew,
                                 double *gpu_sumRref, double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_RecipEnergies, int imageSize);

__global__ void ChangeLambdaMolReciprocalGPU(
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_kx, double *gpu_ky,
    double *gpu_kz, int atomNumber, double *gpu_molCharge, double *gpu_sumRnew,
    double *gpu_sumInew, double *gpu_sumRref, double *gpu_sumIref,
    double *gpu_prefactRef, double *gpu_RecipEnergies, double lambdaCoef,
    int imageSize);

__global__ void SwapReciprocalGPU(
    const double *gpu_x, const double *gpu_y, const double *gpu_z,
    const double *gpu_kx, const double *gpu_ky, const double *gpu_kz,
    int atomNumber, const double *gpu_molCharge, double *gpu_sumRnew,
    double *gpu_sumInew, const double *gpu_sumRref, const double *gpu_sumIref,
    const double *gpu_prefactRef, double *gpu_RecipEnergies, int imageSize);

__global__ void MolExchangeReciprocalGPU(
    int imageSize, double *gpu_kx, double *gpu_ky, double *gpu_kz,
    double *gpu_prefactRef, double *gpu_sumRref, double *gpu_sumIref,
    double *gpu_sumRnew, double *gpu_sumInew, double *gpu_molCharge,
    int numChargedParticles, double *gpu_MolX, double *gpu_MolY,
    double *gpu_MolZ, double *gpu_RecipEnergies, bool first_call);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL_H*/
