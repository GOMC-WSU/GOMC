/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef CALCULATE_ENERGY_CUDA_KERNEL_H
#define CALCULATE_ENERGY_CUDA_KERNEL_H

#ifdef GOMC_CUDA
#include "BoxDimensions.h"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>

void CallBoxInterGPU(VariablesCUDA *vars, const std::vector<int> &cellVector,
                     const std::vector<int> &cellStartIndex,
                     XYZArray const &coords, BoxDimensions const &boxAxes,
                     bool electrostatic, double &REn, double &LJEn,
                     bool sc_coul, double sc_sigma_6, double sc_alpha,
                     uint sc_power, uint const box);

void CallCalculateTorqueGPU(VariablesCUDA *vars,
                            std::vector<int> &moleculeIndex,
                            XYZArray const &coordinates, XYZArray const &com,
                            XYZArray const &atomForce,
                            XYZArray const &atomForceRec, XYZArray &molTorque,
                            const uint box, BoxDimensions const &boxAxes);

__global__ void CalculateTorqueGPU(
    const int *__restrict__ gpu_startAtomIndex, // atom indices
    const int *__restrict__ gpu_moleculeIndex,  // molecule indices

    const double *__restrict__ gpu_cx, // x coordinates
    const double *__restrict__ gpu_cy, // y coordinates
    const double *__restrict__ gpu_cz, // z coordinates

    const double *__restrict__ gpu_comx, // x center of mass
    const double *__restrict__ gpu_comy, // y center of mass
    const double *__restrict__ gpu_comz, // z center of mass

    const double *__restrict__ gpu_atomforcex, // x atom force
    const double *__restrict__ gpu_atomforcey, // y atom force
    const double *__restrict__ gpu_atomforcez, // z atom force

    const double *__restrict__ gpu_atomforcerecx, // x atom force reciprocal
    const double *__restrict__ gpu_atomforcerecy, // y atom force reciprocal
    const double *__restrict__ gpu_atomforcerecz, // z atom force reciprocal

    double *__restrict__ gpu_moltorquex, // x molecule torque
    double *__restrict__ gpu_moltorquey, // y molecule torque
    double *__restrict__ gpu_moltorquez, // z molecule torque

    const uint box,         // box number
    const uint numMolecules, // number of molecules, if index is above this
                            // number we should just return
    const double3 axis,     // lengths in each dimension
    const double3 halfAx    // half of axis
);


__device__ inline double3 CrossProductGPU(const double3 &a, const double3 &b);

__global__ void
BoxInterGPU(int *gpu_cellStartIndex, int *gpu_cellVector, int *gpu_neighborList,
            int numberOfCells, double *gpu_x, double *gpu_y, double *gpu_z,
            double3 axis, double3 halfAx, bool electrostatic,
            double *gpu_particleCharge, int *gpu_particleKind,
            int *gpu_particleMol, double *gpu_REn, double *gpu_LJEn,
            double *gpu_sigmaSq, double *gpu_epsilon_Cn, double *gpu_n,
            int *gpu_VDW_Kind, int *gpu_isMartini, int *gpu_count,
            double *gpu_rCut, double *gpu_rCutSq, double *gpu_rCutCoulomb,
            double *gpu_rCutCoulombSq, double *gpu_rCutLow, double *gpu_rOn,
            double *gpu_alpha, int *gpu_ewald, double *gpu_diElectric_1,
            int *gpu_nonOrth, double *gpu_cell_x, double *gpu_cell_y,
            double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
            double *gpu_Invcell_z, bool sc_coul, double sc_sigma_6,
            double sc_alpha, uint sc_power, double *gpu_rMin,
            double *gpu_rMaxSq, double *gpu_expConst, int *gpu_molIndex,
            double *gpu_lambdaVDW, double *gpu_lambdaCoulomb,
            bool *gpu_isFraction, int box);

__device__ double CalcCoulombGPU(
    double distSq, int kind1, int kind2, double qi_qj_fact, double gpu_rCutLow,
    int gpu_ewald, int gpu_VDW_Kind, double gpu_alpha, double gpu_rCutCoulomb,
    double gpu_rCutCoulombSq, int gpu_isMartini, double gpu_diElectric_1,
    double gpu_lambdaCoulomb, bool sc_coul, double sc_sigma_6, double sc_alpha,
    uint sc_power, double *gpu_sigmaSq, int gpu_count);
__device__ double CalcCoulombVirGPU(double distSq, double qi_qj,
                                    double gpu_rCutCoulomb, double gpu_alpha,
                                    int gpu_VDW_Kind, int gpu_ewald,
                                    double gpu_diElectric_1, int gpu_isMartini);
__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut,
                            double gpu_rCutSq, double gpu_rOn, int gpu_count,
                            double gpu_lambdaVDW, double sc_sigma_6,
                            double sc_alpha, uint sc_power, double *gpu_rMin,
                            double *gpu_rMaxSq, double *gpu_expConst);

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
__device__ double CalcCoulombSwitchMartiniGPU(
    double distSq, int index, double qi_qj_fact, int gpu_ewald,
    double gpu_alpha, double gpu_rCut, double gpu_rCutSq,
    double gpu_diElectric_1, double gpu_lambdaCoulomb, bool sc_coul,
    double sc_sigma_6, double sc_alpha, uint sc_power, double *gpu_sigmaSq);
__device__ double CalcCoulombSwitchMartiniGPUNoLambda(
    double distSq, double qi_qj_fact, int gpu_ewald, double gpu_alpha,
    double gpu_rCut, double gpu_rCutSq, double gpu_diElectric_1);
__device__ double CalcCoulombSwitchGPU(double distSq, int index,
                                       double qi_qj_fact, double gpu_alpha,
                                       int gpu_ewald, double gpu_rCut,
                                       double gpu_rCutSq,
                                       double gpu_lambdaCoulomb, bool sc_coul,
                                       double sc_sigma_6, double sc_alpha,
                                       uint sc_power, double *gpu_sigmaSq);
__device__ double CalcCoulombSwitchGPUNoLambda(double distSq, double qi_qj_fact,
                                               int gpu_ewald, double gpu_alpha,
                                               double gpu_rCut,
                                               double gpu_rCutSq);

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
                                 double gpu_rCutSq, double gpu_lambdaVDW,
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
