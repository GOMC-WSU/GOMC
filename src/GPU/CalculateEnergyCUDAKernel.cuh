/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
                     vector<real> particleCharge,
                     vector<int> particleKind,
                     vector<int> particleMol,
                     real &REn,
                     real &LJEn,
                     bool multiParticleEnabled,
                     real *aForcex,
                     real *aForcey,
                     real *aForcez,
                     real *mForcex,
                     real *mForcey,
                     real *mForcez,
                     int atomCount,
                     int molCount,
                     bool reset_force,
                     bool copy_back,
                     uint const box);

__global__ void BoxInterGPU(int *gpu_pair1,
                            int *gpu_pair2,
                            real *gpu_x,
                            real *gpu_y,
                            real *gpu_z,
                            real xAxes,
                            real yAxes,
                            real zAxes,
                            bool electrostatic,
                            real *gpu_particleCharge,
                            int *gpu_particleKind,
                            int *gpu_particleMol,
                            real *gpu_REn,
                            real *gpu_LJEn,
                            int pairSize,
                            real *gpu_sigmaSq,
                            real *gpu_epsilon_Cn,
                            real *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            real *gpu_rCut,
                            real *gpu_rCutCoulomb,
                            real *gpu_rCutLow,
                            real *gpu_rOn,
                            real *gpu_alpha,
                            int *gpu_ewald,
                            real *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            real *gpu_cell_x,
                            real *gpu_cell_y,
                            real *gpu_cell_z,
                            real *gpu_Invcell_x,
                            real *gpu_Invcell_y,
                            real *gpu_Invcell_z,
                            bool multiParticleEnabled,
                            real *gpu_aForcex,
                            real *gpu_aForcey,
                            real *gpu_aForcez,
                            real *gpu_mForcex,
                            real *gpu_mForcey,
                            real *gpu_mForcez,
                            int box);


__device__ real CalcCoulombGPU(real distSq, real qi_qj_fact,
                                 real gpu_rCutLow, int gpu_ewald,
                                 int gpu_VDW_Kind, real gpu_alpha,
                                 real gpu_rCutCoulomb, int gpu_isMartini,
                                 real gpu_diElectric_1);
__device__ real CalcCoulombVirGPU(real distSq, real qi_qj,
                                    real gpu_rCutCoulomb, real gpu_alpha,
                                    int gpu_VDW_Kind, int gpu_ewald,
                                    real gpu_diElectric_1, int gpu_isMartini);
__device__ real CalcEnGPU(real distSq, int kind1, int kind2,
                            real *gpu_sigmaSq, real *gpu_n,
                            real *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, real gpu_rCut, real gpu_rOn,
                            int gpu_count);

//ElectroStatic Calculation
//**************************************************************//
__device__ real CalcCoulombParticleGPU(real distSq, real qi_qj_fact,
                                         real gpu_ewald, real gpu_alpha);
__device__ real CalcCoulombShiftGPU(real distSq, real qi_qj_fact,
                                      int gpu_ewald, real gpu_alpha,
                                      real gpu_rCut);
__device__ real CalcCoulombSwitchMartiniGPU(real distSq, real qi_qj_fact,
                                              int gpu_ewald, real gpu_alpha,
                                              real gpu_rCut,
                                              real gpu_diElectric_1);
__device__ real CalcCoulombSwitchGPU(real distSq, real qi_qj_fact,
                                       real gpu_alpha, int gpu_ewald,
                                       real gpu_rCut);

//VDW Calculation
//*****************************************************************//
__device__ real CalcEnParticleGPU(real distSq, int index,
                                    real *gpu_sigmaSq, real *gpu_n,
                                    real *gpu_epsilon_Cn);
__device__ real CalcEnShiftGPU(real distSq, int index, real *gpu_sigmaSq,
                                 real *gpu_n, real *gpu_epsilon_Cn,
                                 real gpu_rCut);
__device__ real CalcEnSwitchMartiniGPU(real distSq, int index,
    real *gpu_sigmaSq, real *gpu_n,
    real *gpu_epsilon_Cn,
    real gpu_rCut, real gpu_rOn);
__device__ real CalcEnSwitchGPU(real distSq, int index, real *gpu_sigmaSq,
                                  real *gpu_n, real *gpu_epsilon_Cn,
                                  real gpu_rCut, real gpu_rOn);

#endif /*GOMC_CUDA*/
