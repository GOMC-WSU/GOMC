/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"
#include "VariablesCUDA.cuh"

using namespace std;

void CallBoxInterGPU(VariablesCUDA *vars,
                     vector<uint> pair1,
                     vector<uint> pair2,
                     XYZArray const &coords,
                     BoxDimensions const &boxAxes,
                     bool electrostatic,
                     vector<double> particleCharge,
                     vector<int> particleKind,
                     vector<int> particleMol,
                     double &REn,
                     double &LJEn,
                     bool sc_coul,
                     double sc_sigma_6,
                     double sc_alpha,
                     uint sc_power,
                     uint const box);

__global__ void BoxInterGPU(int *gpu_pair1,
                            int *gpu_pair2,
                            double *gpu_x,
                            double *gpu_y,
                            double *gpu_z,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            bool electrostatic,
                            double *gpu_particleCharge,
                            int *gpu_particleKind,
                            int *gpu_particleMol,
                            double *gpu_REn,
                            double *gpu_LJEn,
                            int pairSize,
                            double *gpu_sigmaSq,
                            double *gpu_epsilon_Cn,
                            double *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            double *gpu_rCut,
                            double *gpu_rCutCoulomb,
                            double *gpu_rCutLow,
                            double *gpu_rOn,
                            double *gpu_alpha,
                            int *gpu_ewald,
                            double *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            double *gpu_cell_x,
                            double *gpu_cell_y,
                            double *gpu_cell_z,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            bool sc_coul,
                            double sc_sigma_6,
                            double sc_alpha,
                            uint sc_power,
                            double *gpu_rMin,
                            double *gpu_rMaxSq,
                            double *gpu_expConst,
                            int box);


__device__ double CalcCoulombGPU(double distSq, int kind1, int kind2,
                                 double qi_qj_fact, double gpu_rCutLow,
                                 int gpu_ewald, int gpu_VDW_Kind,
                                 double gpu_alpha, double gpu_rCutCoulomb,
                                 int gpu_isMartini, double gpu_diElectric_1,
                                 bool sc_coul,
                                 double sc_sigma_6, double sc_alpha,
                                 uint sc_power, double gpu_sigmaSq,
                                 int gpu_count);
__device__ double CalcCoulombVirGPU(double distSq, double qi_qj,
                                    double gpu_rCutCoulomb, double gpu_alpha,
                                    int gpu_VDW_Kind, int gpu_ewald,
                                    double gpu_diElectric_1, int gpu_isMartini);
__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut, double gpu_rOn,
                            int gpu_count, double sc_sigma_6, double sc_alpha,
                            uint sc_power, double *gpu_rMin, double *gpu_rMaxSq,
                            double *gpu_expConst);

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact,
    double gpu_ewald, double gpu_alpha, bool sc_coul,
    double sc_sigma_6, double sc_alpha,
    uint sc_power, double gpu_sigmaSq);

__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj,
    double gpu_ewald, double gpu_alpha);

//VDW Calculation
//*****************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index,
                                    double *gpu_sigmaSq, double *gpu_n,
                                    double *gpu_epsilon_Cn,
                                    double sc_sigma_6,
                                    double sc_alpha,
                                    uint sc_power);

#endif /*GOMC_CUDA*/
