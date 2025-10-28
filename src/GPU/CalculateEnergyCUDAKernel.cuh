/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CALCULATE_ENERGY_CUDA_KERNEL_H
#define CALCULATE_ENERGY_CUDA_KERNEL_H

#ifdef GOMC_CUDA
#include "BoxDimensions.h"
#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>

void CallBoxInterGPU(VariablesCUDA *vars, const std::vector<int> &cellVector,
                     const std::vector<int> &cellStartIndex,
                     const std::vector<std::vector<int>> &neighborList,
                     XYZArray const &coords, BoxDimensions const &boxAxes,
                     bool electrostatic,
                     const std::vector<double> &particleCharge,
                     const std::vector<int> &particleKind,
                     const std::vector<int> &particleMol, double &REn,
                     double &LJEn, bool sc_coul, double sc_sigma_6,
                     double sc_alpha, uint sc_power, uint const box);

__global__ void
BoxInterGPU(int *gpu_cellStartIndex, int *gpu_cellVector, int *gpu_neighborList,
            int numberOfCells, double *gpu_x, double *gpu_y, double *gpu_z,
            double3 axis, double3 halfAx, bool electrostatic,
            double *gpu_particleCharge, int *gpu_particleKind,
            int *gpu_particleMol, double *gpu_REn, double *gpu_LJEn,
            double *gpu_sigmaSq, double *gpu_epsilon_Cn, double *gpu_n,
            int *gpu_VDW_Kind, int *gpu_isMartini, int *gpu_count,
            double *gpu_rCut, double *gpu_rCutCoulomb, double *gpu_rCutLow,
            double *gpu_rOn, double *gpu_alpha, int *gpu_ewald,
            double *gpu_diElectric_1, int *gpu_nonOrth, double *gpu_cell_x,
            double *gpu_cell_y, double *gpu_cell_z, double *gpu_Invcell_x,
            double *gpu_Invcell_y, double *gpu_Invcell_z, bool sc_coul,
            double sc_sigma_6, double sc_alpha, uint sc_power, double *gpu_rMin,
            double *gpu_rMaxSq, double *gpu_expConst, int *gpu_molIndex,
            double *gpu_lambdaVDW, double *gpu_lambdaCoulomb,
            bool *gpu_isFraction, int box);

__device__ double
CalcCoulombGPU(double distSq, int kind1, int kind2, double qi_qj_fact,
               double gpu_rCutLow, int gpu_ewald, int gpu_VDW_Kind,
               double gpu_alpha, double gpu_rCutCoulomb, int gpu_isMartini,
               double gpu_diElectric_1, double gpu_lambdaCoulomb, bool sc_coul,
               double sc_sigma_6, double sc_alpha, uint sc_power,
               double *gpu_sigmaSq, int gpu_count);
__device__ double CalcCoulombVirGPU(double distSq, double qi_qj,
                                    double gpu_rCutCoulomb, double gpu_alpha,
                                    int gpu_VDW_Kind, int gpu_ewald,
                                    double gpu_diElectric_1, int gpu_isMartini);
__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut, double gpu_rOn,
                            int gpu_count, double gpu_lambdaVDW,
                            double sc_sigma_6, double sc_alpha, uint sc_power,
                            double *gpu_rMin, double *gpu_rMaxSq,
                            double *gpu_expConst);

// ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, int index,
                                         double qi_qj_fact, int gpu_ewald,
                                         double gpu_alpha,
                                         double gpu_lambdaCoulomb, bool sc_coul,
                                         double sc_sigma_6, double sc_alpha,
                                         uint sc_power, double *gpu_sigmaSq);
__device__ double CalcCoulombParticleGPUNoLambda(double distSq,
                                                 double qi_qj_fact,
                                                 int gpu_ewald,
                                                 double gpu_alpha);
__device__ double CalcCoulombShiftGPU(double distSq, int index,
                                      double qi_qj_fact, int gpu_ewald,
                                      double gpu_alpha, double gpu_rCut,
                                      double gpu_lambdaCoulomb, bool sc_coul,
                                      double sc_sigma_6, double sc_alpha,
                                      uint sc_power, double *gpu_sigmaSq);
__device__ double CalcCoulombShiftGPUNoLambda(double distSq, double qi_qj_fact,
                                              int gpu_ewald, double gpu_alpha,
                                              double gpu_rCut);
__device__ double CalcCoulombExp6GPU(double distSq, int index,
                                     double qi_qj_fact, int gpu_ewald,
                                     double gpu_alpha, double gpu_lambdaCoulomb,
                                     bool sc_coul, double sc_sigma_6,
                                     double sc_alpha, uint sc_power,
                                     double *gpu_sigmaSq);
__device__ double CalcCoulombExp6GPUNoLambda(double distSq, double qi_qj_fact,
                                             int gpu_ewald, double gpu_alpha);
__device__ double
CalcCoulombSwitchMartiniGPU(double distSq, int index, double qi_qj_fact,
                            int gpu_ewald, double gpu_alpha, double gpu_rCut,
                            double gpu_diElectric_1, double gpu_lambdaCoulomb,
                            bool sc_coul, double sc_sigma_6, double sc_alpha,
                            uint sc_power, double *gpu_sigmaSq);
__device__ double
CalcCoulombSwitchMartiniGPUNoLambda(double distSq, double qi_qj_fact,
                                    int gpu_ewald, double gpu_alpha,
                                    double gpu_rCut, double gpu_diElectric_1);
__device__ double CalcCoulombSwitchGPU(double distSq, int index,
                                       double qi_qj_fact, double gpu_alpha,
                                       int gpu_ewald, double gpu_rCut,
                                       double gpu_lambdaCoulomb, bool sc_coul,
                                       double sc_sigma_6, double sc_alpha,
                                       uint sc_power, double *gpu_sigmaSq);
__device__ double CalcCoulombSwitchGPUNoLambda(double distSq, double qi_qj_fact,
                                               int gpu_ewald, double gpu_alpha,
                                               double gpu_rCut);

// VDW Calculation
//*****************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index,
                                    double *gpu_sigmaSq, double *gpu_n,
                                    double *gpu_epsilon_Cn,
                                    double gpu_lambdaVDW, double sc_sigma_6,
                                    double sc_alpha, uint sc_power);
__device__ double CalcEnParticleGPUNoLambda(double distSq, int index,
                                            double *gpu_sigmaSq, double *gpu_n,
                                            double *gpu_epsilon_Cn);
__device__ double CalcEnShiftGPU(double distSq, int index, double *gpu_sigmaSq,
                                 double *gpu_n, double *gpu_epsilon_Cn,
                                 double gpu_rCut, double gpu_lambdaVDW,
                                 double sc_sigma_6, double sc_alpha,
                                 uint sc_power);
__device__ double CalcEnShiftGPUNoLambda(double distSq, int index,
                                         double *gpu_sigmaSq, double *gpu_n,
                                         double *gpu_epsilon_Cn,
                                         double gpu_rCut);
__device__ double CalcEnExp6GPU(double distSq, int index, double *gpu_sigmaSq,
                                double *gpu_n, double gpu_lambdaVDW,
                                double sc_sigma_6, double sc_alpha,
                                uint sc_power, double *gpu_rMin,
                                double *gpu_rMaxSq, double *gpu_expConst);
__device__ double CalcEnExp6GPUNoLambda(double distSq, int index, double *gpu_n,
                                        double *gpu_rMin, double *gpu_expConst);
__device__ double
CalcEnSwitchMartiniGPU(double distSq, int index, double *gpu_sigmaSq,
                       double *gpu_n, double *gpu_epsilon_Cn, double gpu_rCut,
                       double gpu_rOn, double gpu_lambdaVDW, double sc_sigma_6,
                       double sc_alpha, uint sc_power);
__device__ double
CalcEnSwitchMartiniGPUNoLambda(double distSq, int index, double *gpu_sigmaSq,
                               double *gpu_n, double *gpu_epsilon_Cn,
                               double gpu_rCut, double gpu_rOn);
__device__ double CalcEnSwitchGPU(double distSq, int index, double *gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn,
                                  double gpu_rCut, double gpu_rOn,
                                  double gpu_lambdaVDW, double sc_sigma_6,
                                  double sc_alpha, uint sc_power);
__device__ double CalcEnSwitchGPUNoLambda(double distSq, int index,
                                          double *gpu_sigmaSq, double *gpu_n,
                                          double *gpu_epsilon_Cn,
                                          double gpu_rCut, double gpu_rOn);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_ENERGY_CUDA_KERNEL_H*/
