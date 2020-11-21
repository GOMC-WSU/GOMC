/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA
#include <vector>
#include "Random123/philox.h"
typedef r123::Philox4x32 RNG;

#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include "XYZArray.h"

void CallTranslateParticlesGPU(VariablesCUDA *vars,
                               std::vector<int> &isParticleInvolved,
                               double t_max,
                               double *mForcex,
                               double *mForcey,
                               double *mForcez,
                               unsigned int step,
                               unsigned int seed,
                               std::vector<int> particleMol,
                               int atomCount,
                               int molCount,
                               double xAxes,
                               double yAxes,
                               double zAxes,
                               XYZArray &newMolPos,
                               XYZArray &newCOMs,
                               double lambdaBETA,
                               XYZArray &t_k,
                               XYZArray &molForceRecRef);

void CallRotateParticlesGPU(VariablesCUDA *vars,
                            std::vector<int> &isParticleInvolved,
                            double r_max,
                            double *mTorquex,
                            double *mTorquey,
                            double *mTorquez,
                            unsigned int step,
                            unsigned int seed,
                            std::vector<int> particleMol,
                            int atomCount,
                            int molCount,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            XYZArray &newMolPos,
                            XYZArray &newCOMs,
                            double lambdaBETA,
                            XYZArray &r_k);

__global__ void TranslateParticlesKernel(unsigned int numberOfMolecules,
    double t_max,
    double *molForcex,
    double *molForcey,
    double *molForcez,
    unsigned int step,
    unsigned int seed,
    double *gpu_x,
    double *gpu_y,
    double *gpu_z,
    int *gpu_particleMol,
    int atomCount,
    double xAxes,
    double yAxes,
    double zAxes,
    double *gpu_comx,
    double *gpu_comy,
    double *gpu_comz,
    double lambdaBETA,
    double *gpu_t_k_x,
    double *gpu_t_k_y,
    double *gpu_t_k_z,
    int *gpu_isParticleInvolved,
    double *gpu_mForceRecx,
    double *gpu_mForceRecy,
    double *gpu_mForceRecz);

__global__ void RotateParticlesKernel(unsigned int numberOfMolecules,
                                      double r_max,
                                      double *molTorquex,
                                      double *molTorquey,
                                      double *molTorquez,
                                      unsigned int step,
                                      unsigned int seed,
                                      double *gpu_x,
                                      double *gpu_y,
                                      double *gpu_z,
                                      int *gpu_particleMol,
                                      int atomCount,
                                      double xAxes,
                                      double yAxes,
                                      double zAxes,
                                      double *gpu_comx,
                                      double *gpu_comy,
                                      double *gpu_comz,
                                      double lambdaBETA,
                                      double *gpu_r_k_x,
                                      double *gpu_r_k_y,
                                      double *gpu_r_k_z,
                                      int *gpu_isParticleInvolved);
#endif
