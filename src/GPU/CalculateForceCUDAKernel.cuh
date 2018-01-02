/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
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

using namespace std;

void CallBoxInterForceGPU(VariablesCUDA *vars,
                          vector<uint> &pair1,
                          vector<uint> &pair2,
                          XYZArray const &currentCoords,
                          XYZArray const &currentCOM,
                          BoxDimensions const& boxAxes,
                          bool electrostatic,
                          vector<double> &particleCharge,
                          vector<int> &particleKind,
                          vector<int> &particleMol,
                          double &rT11,
                          double &rT12,
                          double &rT13,
                          double &rT22,
                          double &rT23,
                          double &rT33,
                          double &vT11,
                          double &vT12,
                          double &vT13,
                          double &vT22,
                          double &vT23,
                          double &vT33,
                          uint const box);

void CallForceReciprocalGPU(VariablesCUDA *vars,
                            XYZArray const &currentCoords,
                            XYZArray const &currentCOMDiff,
                            vector<double> &particleCharge,
                            double &rT11,
                            double &rT12,
                            double &rT13,
                            double &rT22,
                            double &rT23,
                            double &rT33,
                            uint imageSize,
                            double constVal,
                            uint box);

__global__ void BoxInterForceGPU(int *gpu_pair1,
                                 int *gpu_pair2,
                                 double *gpu_x,
                                 double *gpu_y,
                                 double *gpu_z,
                                 double *gpu_comx,
                                 double *gpu_comy,
                                 double *gpu_comz,
                                 double xAxes,
                                 double yAxes,
                                 double zAxes,
                                 bool electrostatic,
                                 double *gpu_particleCharge,
                                 int *gpu_particleKind,
                                 int *gpu_particleMol,
                                 double *gpu_rT11,
                                 double *gpu_rT12,
                                 double *gpu_rT13,
                                 double *gpu_rT22,
                                 double *gpu_rT23,
                                 double *gpu_rT33,
                                 double *gpu_vT11,
                                 double *gpu_vT12,
                                 double *gpu_vT13,
                                 double *gpu_vT22,
                                 double *gpu_vT23,
                                 double *gpu_vT33,
                                 int pairSize,
                                 double *gpu_sigmaSq,
                                 double *gpu_epsilon_Cn,
                                 double *gpu_n,
                                 int *gpu_VDW_Kind,
                                 int *gpu_isMartini,
                                 int *gpu_count,
                                 double *gpu_rCut,
                                 double *gpu_rCutLow,
                                 double *gpu_rOn,
                                 double *gpu_alpha,
                                 int *gpu_ewald,
                                 double *gpu_diElectric_1);

__global__ void ForceReciprocalGPU(double *gpu_x,
                                   double *gpu_y,
                                   double *gpu_z,
                                   double *gpu_comDx,
                                   double *gpu_comDy,
                                   double *gpu_comDz,
                                   double *gpu_kxRef,
                                   double *gpu_kyRef,
                                   double *gpu_kzRef,
                                   double *gpu_prefactRef,
                                   double *gpu_hsqrRef,
                                   double *gpu_sumRref,
                                   double *gpu_sumIref,
                                   double *gpu_particleCharge,
                                   double *gpu_rT11,
                                   double *gpu_rT12,
                                   double *gpu_rT13,
                                   double *gpu_rT22,
                                   double *gpu_rT23,
                                   double *gpu_rT33,
                                   double constVal,
                                   uint imageSize,
                                   uint atomNumber);

__device__ double CalcCoulombForceGPU(double distSq, double qi_qj,
                                      int gpu_VDW_Kind,
                                      int gpu_ewald,
                                      int gpu_isMartini,
                                      double gpu_alpha,
                                      double gpu_rCut,
                                      double gpu_diElectric_1);
__device__ double CalcEnForceGPU(double distSq, int kind1, int kind2,
                                 double *gpu_sigmaSq,
                                 double *gpu_n,
                                 double *gpu_epsilon_Cn,
                                 double gpu_rCut,
                                 double gpu_rOn,
                                 int gpu_isMartini,
                                 int gpu_VDW_Kind,
                                 int gpu_count);

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj,
    double gpu_alpha);
__device__ double CalcCoulombVirShiftGPU(double distSq, double qi_qj,
    int gpu_ewald, double gpu_alpha);
__device__ double CalcCoulombVirSwitchMartiniGPU(double distSq, double qi_qj,
    int gpu_ewald,
    double gpu_alpha,
    double gpu_rCut,
    double gpu_diElectric_1);
__device__ double CalcCoulombVirSwitchGPU(double distSq, double qi_qj,
    int gpu_ewald, double gpu_alpha,
    double gpu_rCut);

//VDW Calculation
//*****************************************************************//
__device__ double CalcVirParticleGPU(double distSq, int index,
                                     double *gpu_sigmaSq, double *gpu_n,
                                     double *gpu_epsilon_Cn);
__device__ double CalcVirShiftGPU(double distSq, int index,
                                  double *gpu_sigmaSq, double *gpu_n,
                                  double *gpu_epsilon_Cn);
__device__ double CalcVirSwitchMartiniGPU(double distSq, int index,
    double *gpu_sigmaSq, double *gpu_n,
    double *gpu_epsilon_Cn,
    double gpu_rCut, double rOn);
__device__ double CalcVirSwitchGPU(double distSq, int index,
                                   double *gpu_sigmaSq, double *gpu_epsilon_Cn,
                                   double *gpu_n, double gpu_rCut,
                                   double gpu_rOn);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_FORCE_CUDA_KERNEL*/
