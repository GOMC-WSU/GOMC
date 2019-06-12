/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CALCULATE_FORCE_CUDA_KERNEL
#define CALCULATE_FORCE_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"
#include "VariablesCUDA.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"

using namespace std;

void CallBoxInterForceGPU(VariablesCUDA *vars,
                          vector<uint> &pair1,
                          vector<uint> &pair2,
                          XYZArray const &currentCoords,
                          XYZArray const &currentCOM,
                          BoxDimensions const& boxAxes,
                          bool electrostatic,
                          vector<real> &particleCharge,
                          vector<int> &particleKind,
                          vector<int> &particleMol,
                          real &rT11,
                          real &rT12,
                          real &rT13,
                          real &rT22,
                          real &rT23,
                          real &rT33,
                          real &vT11,
                          real &vT12,
                          real &vT13,
                          real &vT22,
                          real &vT23,
                          real &vT33,
                          uint const box);

void CallForceReciprocalGPU(VariablesCUDA *vars,
                            XYZArray const &currentCoords,
                            XYZArray const &currentCOMDiff,
                            vector<real> &particleCharge,
                            real &rT11,
                            real &rT12,
                            real &rT13,
                            real &rT22,
                            real &rT23,
                            real &rT33,
                            uint imageSize,
                            real constVal,
                            uint box);

__global__ void BoxInterForceGPU(int *gpu_pair1,
                                 int *gpu_pair2,
                                 real *gpu_x,
                                 real *gpu_y,
                                 real *gpu_z,
                                 real *gpu_comx,
                                 real *gpu_comy,
                                 real *gpu_comz,
                                 real xAxes,
                                 real yAxes,
                                 real zAxes,
                                 bool electrostatic,
                                 real *gpu_particleCharge,
                                 int *gpu_particleKind,
                                 int *gpu_particleMol,
                                 real *gpu_rT11,
                                 real *gpu_rT12,
                                 real *gpu_rT13,
                                 real *gpu_rT22,
                                 real *gpu_rT23,
                                 real *gpu_rT33,
                                 real *gpu_vT11,
                                 real *gpu_vT12,
                                 real *gpu_vT13,
                                 real *gpu_vT22,
                                 real *gpu_vT23,
                                 real *gpu_vT33,
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
                                 real *gpu_cell_x,
                                 real *gpu_cell_y,
                                 real *gpu_cell_z,
                                 real *gpu_Invcell_x,
                                 real *gpu_Invcell_y,
                                 real *gpu_Invcell_z,
                                 int *gpu_nonOrth,
                                 int box);

__global__ void ForceReciprocalGPU(real *gpu_x,
                                   real *gpu_y,
                                   real *gpu_z,
                                   real *gpu_comDx,
                                   real *gpu_comDy,
                                   real *gpu_comDz,
                                   real *gpu_kxRef,
                                   real *gpu_kyRef,
                                   real *gpu_kzRef,
                                   real *gpu_prefactRef,
                                   real *gpu_hsqrRef,
                                   real *gpu_sumRref,
                                   real *gpu_sumIref,
                                   real *gpu_particleCharge,
                                   real *gpu_rT11,
                                   real *gpu_rT12,
                                   real *gpu_rT13,
                                   real *gpu_rT22,
                                   real *gpu_rT23,
                                   real *gpu_rT33,
                                   real constVal,
                                   uint imageSize,
                                   uint atomNumber);

__device__ real CalcCoulombForceGPU(real distSq, real qi_qj,
                                      int gpu_VDW_Kind,
                                      int gpu_ewald,
                                      int gpu_isMartini,
                                      real gpu_alpha,
                                      real gpu_rCutCoulomb,
                                      real gpu_diElectric_1);
__device__ real CalcEnForceGPU(real distSq, int kind1, int kind2,
                                 real *gpu_sigmaSq,
                                 real *gpu_n,
                                 real *gpu_epsilon_Cn,
                                 real gpu_rCut,
                                 real gpu_rOn,
                                 int gpu_isMartini,
                                 int gpu_VDW_Kind,
                                 int gpu_count);

//ElectroStatic Calculation
//**************************************************************//
__device__ real CalcCoulombVirParticleGPU(real distSq, real qi_qj,
    real gpu_alpha);
__device__ real CalcCoulombVirShiftGPU(real distSq, real qi_qj,
    int gpu_ewald, real gpu_alpha);
__device__ real CalcCoulombVirSwitchMartiniGPU(real distSq, real qi_qj,
    int gpu_ewald,
    real gpu_alpha,
    real gpu_rCut,
    real gpu_diElectric_1);
__device__ real CalcCoulombVirSwitchGPU(real distSq, real qi_qj,
    int gpu_ewald, real gpu_alpha,
    real gpu_rCut);

//VDW Calculation
//*****************************************************************//
__device__ real CalcVirParticleGPU(real distSq, int index,
                                     real *gpu_sigmaSq, real *gpu_n,
                                     real *gpu_epsilon_Cn);
__device__ real CalcVirShiftGPU(real distSq, int index,
                                  real *gpu_sigmaSq, real *gpu_n,
                                  real *gpu_epsilon_Cn);
__device__ real CalcVirSwitchMartiniGPU(real distSq, int index,
    real *gpu_sigmaSq, real *gpu_n,
    real *gpu_epsilon_Cn,
    real gpu_rCut, real rOn);
__device__ real CalcVirSwitchGPU(real distSq, int index,
                                   real *gpu_sigmaSq, real *gpu_epsilon_Cn,
                                   real *gpu_n, real gpu_rCut,
                                   real gpu_rOn);


// Have to move the implementation for some functions here 
// since CUDA doesn't allow __global__ to call __device__
// from different files
// Wanted to call CalcCoulombForceGPU() from CalculateEnergyCUDAKernel.cu file
__device__ inline real CalcCoulombForceGPU(real distSq, real qi_qj,
  int gpu_VDW_Kind, int gpu_ewald,
  int gpu_isMartini, real gpu_alpha,
  real gpu_rCutCoulomb,
  real gpu_diElectric_1)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombVirParticleGPU(distSq, qi_qj, gpu_alpha);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcCoulombVirShiftGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                          gpu_rCutCoulomb, gpu_diElectric_1);
  } else
    return CalcCoulombVirSwitchGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                   gpu_rCutCoulomb);
}


#endif /*GOMC_CUDA*/
#endif /*CALCULATE_FORCE_CUDA_KERNEL*/
