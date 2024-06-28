/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CALCULATE_EWALD_CUDA_KERNEL
#define CALCULATE_EWALD_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include <vector>
#include "XYZArray.h"

void CallBoxForceReciprocalGPU(VariablesCUDA *vars,
                               XYZArray &atomForceRec,
                               XYZArray &molForceRec,
                               const std::vector<double> &particleCharge,
                               const std::vector<int> &particleMol,
                               const std::vector<bool> &particleHasNoCharge,
                               const bool *particleUsed,
                               const std::vector<int> &startMol,
                               const std::vector<int> &lengthMol,
                               double alpha,
                               double alphaSq,
                               double constValue,
                               uint imageSize,
                               XYZArray const &molCoords,
                               BoxDimensions const &boxAxes,
                               int box);

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars,
                               XYZArray const &coords,
                               double const *kx,
                               double const *ky,
                               double const *kz,
                               const std::vector<double> &particleCharge,
                               uint imageSize,
                               double *sumRnew,
                               double *sumInew,
                               double *prefact,
                               double *hsqr,
                               double &energyRecip,
                               uint box);

void CallBoxReciprocalSumsGPU(VariablesCUDA *vars,
                              XYZArray const &coords,
                              const std::vector<double> &particleCharge,
                              uint imageSize,
                              double *sumRnew,
                              double *sumInew,
                              double &energyRecip,
                              uint box);

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          const std::vector<double> &particleCharge,
                          uint imageSize,
                          double *sumRnew,
                          double *sumInew,
                          double &energyRecip,
                          uint box);

//Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &particleCharge,
                                      uint imageSize,
                                      double *sumRnew,
                                      double *sumInew,
                                      double &energyRecip,
                                      const double lambdaCoef,
                                      uint box);

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           const std::vector<double> &particleCharge,
                           uint imageSize,
                           double *sumRnew,
                           double *sumInew,
                           double *sumRref,
                           double *sumIref,
                           const bool insert,
                           double &energyRecip,
                           uint box);

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars,
                                  uint imageSize,
                                  double *sumRnew,
                                  double *sumInew,
                                  uint box, 
                                  std::vector<double> chargeBoxNew,
                                  std::vector<double> chargeBoxOld,
                                  uint lengthNew, uint lengthOld,
                                  double &energyRecipNew,
                                  XYZArray newMolCoords,
                                  XYZArray oldMolCoords);


__global__ void BoxForceReciprocalGPU(double *gpu_aForceRecx,
                                      double *gpu_aForceRecy,
                                      double *gpu_aForceRecz,
                                      double *gpu_mForceRecx,
                                      double *gpu_mForceRecy,
                                      double *gpu_mForceRecz,
                                      double *gpu_particleCharge,
                                      int *gpu_particleMol,
                                      bool *gpu_particleHasNoCharge,
                                       bool *gpu_particleUsed,
                                     int *gpu_startMol,
                                      int *gpu_lengthMol,
                                      double alpha,
                                      double alphaSq,
                                      double constValue,
                                      int imageSize,
                                      double *gpu_kx,
                                      double *gpu_ky,
                                      double *gpu_kz,
                                      double *gpu_x,
                                      double *gpu_y,
                                      double *gpu_z,
                                      double *gpu_prefact,
                                      double *gpu_sumRnew,
                                      double *gpu_sumInew,
                                      bool *gpu_isFraction,
                                      int *gpu_molIndex,
                                      double *gpu_lambdaCoulomb,
                                      double *gpu_cell_x,
                                      double *gpu_cell_y,
                                      double *gpu_cell_z,
                                      double *gpu_Invcell_x,
                                      double *gpu_Invcell_y,
                                      double *gpu_Invcell_z,
                                      int *gpu_nonOrth,
                                      double axx,
                                      double axy,
                                      double axz,
                                      int box,
                                      int atomCount);

__global__ void BoxReciprocalSumsGPU(double * gpu_x,
                                     double * gpu_y,
                                     double * gpu_z,
                                     double * gpu_kx,
                                     double * gpu_ky,
                                     double * gpu_kz,
                                     int atomNumber,
                                     double * gpu_particleCharge,
                                     double * gpu_sumRnew,
                                     double * gpu_sumInew,
                                     int imageSize);

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber,
                                 double *gpu_particleCharge,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_sumRref,
                                 double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_RecipEnergies,
                                 int imageSize);

__global__ void ChangeLambdaMolReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                            double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                            int atomNumber,
                                            double *gpu_particleCharge,
                                            double *gpu_sumRnew,
                                            double *gpu_sumInew,
                                            double *gpu_sumRref,
                                            double *gpu_sumIref,
                                            double *gpu_prefactRef,
                                            double *gpu_RecipEnergies,
                                            double lambdaCoef,
                                            int imageSize);

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                  double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                  int atomNumber,
                                  double *gpu_particleCharge,
                                  double *gpu_sumRnew,
                                  double *gpu_sumInew,
                                  double *gpu_sumRref,
                                  double *gpu_sumIref,
                                  double *gpu_prefactRef,
                                  const bool insert,
                                  double *gpu_RecipEnergies,
                                  int imageSize);

__global__ void MolExchangeReciprocalGPU(
                                  int imageSize,
                                  double *gpu_kx,
                                  double *gpu_ky,
                                  double *gpu_kz,
                                  double *gpu_sumRnew,
                                  double *gpu_sumInew,
                                  double *gpu_chargeBoxNew,
                                  double *gpu_chargeBoxOld,
                                  uint lengthNew,
                                  uint lengthOld,
                                  double *gpu_newMolX,
                                  double *gpu_newMolY,
                                  double *gpu_newMolZ,
                                  double *gpu_oldMolX,
                                  double *gpu_oldMolY,
                                  double *gpu_oldMolZ);

__global__ void NewSwapReciprocalGPU(VariablesCUDA *vars,
                                  int atomNumber,
                                  uint box,
                                  double *gpu_particleCharge,
                                  const bool insert,
                                  double *gpu_RecipEnergies,
                                  int imageSize);

__global__ void BoxReciprocalGPU(double *gpu_prefact,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_RecipEnergies,
                                 int imageSize);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL*/
