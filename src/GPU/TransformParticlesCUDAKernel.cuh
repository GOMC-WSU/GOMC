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
typedef r123::Philox4x64 RNG;
static const double RAND_INTERVAL_GPU = 1.0/static_cast<double>(ULONG_MAX);

#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include "math.h"

void CallTranslateParticlesGPU(VariablesCUDA *vars,
                               const std::vector<int8_t> &isMoleculeInvolved,
                               int box,
                               double t_max,
                               double *mForcex,
                               double *mForcey,
                               double *mForcez,
                               ulong step,
                               ulong seed,
                               const std::vector<int> &particleMol,
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
                            const std::vector<int8_t> &isMoleculeInvolved,
                            int box,
                            double r_max,
                            double *mTorquex,
                            double *mTorquey,
                            double *mTorquez,
                            ulong step,
                            ulong seed,
                            const std::vector<int> &particleMol,
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
    ulong step,
    ulong seed,
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
    double *gpu_cell_x,
    double *gpu_cell_y,
    double *gpu_cell_z,
    double *gpu_Invcell_x,
    double *gpu_Invcell_y,
    double *gpu_Invcell_z,
    int *gpu_nonOrth,
    double lambdaBETA,
    double *gpu_t_k_x,
    double *gpu_t_k_y,
    double *gpu_t_k_z,
    int8_t *gpu_isMoleculeInvolved,
    double *gpu_mForceRecx,
    double *gpu_mForceRecy,
    double *gpu_mForceRecz);

__global__ void RotateParticlesKernel(unsigned int numberOfMolecules,
                                      double r_max,
                                      double *molTorquex,
                                      double *molTorquey,
                                      double *molTorquez,
                                      ulong step,
                                      ulong seed,
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
                                      double *gpu_cell_x,
                                      double *gpu_cell_y,
                                      double *gpu_cell_z,
                                      double *gpu_Invcell_x,
                                      double *gpu_Invcell_y,
                                      double *gpu_Invcell_z,
                                      int *gpu_nonOrth,
                                      double lambdaBETA,
                                      double *gpu_r_k_x,
                                      double *gpu_r_k_y,
                                      double *gpu_r_k_z,
                                      int8_t *gpu_isMoleculeInvolved);

// Brownian Motion multiparticle
void BrownianMotionRotateParticlesGPU(
  VariablesCUDA *vars,
  const std::vector<unsigned int> &moleculeInvolved,
  XYZArray &mTorque,
  XYZArray &newMolPos,
  XYZArray &newCOMs,
  XYZArray &r_k,
  const XYZ &boxAxes,
  const double BETA,
  const double r_max,
  ulong step,
  ulong seed,
  const int box,
  const bool isOrthogonal,
  int *kill);


void BrownianMotionTranslateParticlesGPU(
  VariablesCUDA *vars,
  const std::vector<unsigned int> &moleculeInvolved,
  XYZArray &mForce,
  XYZArray &mForceRec,
  XYZArray &newMolPos,
  XYZArray &newCOMs,
  XYZArray &t_k,
  const XYZ &boxAxes,
  const double BETA,
  const double t_max,
  ulong step,
  ulong seed,
  const int box,
  const bool isOrthogonal,
  int *kill);


template<const bool isOrthogonal>
__global__ void BrownianMotionRotateKernel(
  int *startAtomIdx,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  double *molTorquex,
  double *molTorquey,
  double *molTorquez,
  double *gpu_comx,
  double *gpu_comy,
  double *gpu_comz,
  double *gpu_r_k_x,
  double *gpu_r_k_y,
  double *gpu_r_k_z,
  int *moleculeInvolved,
  double *gpu_cell_x,
  double *gpu_cell_y,
  double *gpu_cell_z,
  double *gpu_Invcell_x,
  double *gpu_Invcell_y,
  double *gpu_Invcell_z,
  double3 axis,
  double3 halfAx,
  int atomCount,
  double r_max,
  ulong step,
  ulong seed,
  double BETA,
  int *kill);


template<const bool isOrthogonal>
__global__ void BrownianMotionTranslateKernel(
  int *startAtomIdx,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  double *molForcex,
  double *molForcey,
  double *molForcez,
  double *molForceRecx,
  double *molForceRecy,
  double *molForceRecz,
  double *gpu_comx,
  double *gpu_comy,
  double *gpu_comz,
  double *gpu_t_k_x,
  double *gpu_t_k_y,
  double *gpu_t_k_z,
  int *moleculeInvolved,
  double *gpu_cell_x,
  double *gpu_cell_y,
  double *gpu_cell_z,
  double *gpu_Invcell_x,
  double *gpu_Invcell_y,
  double *gpu_Invcell_z,
  double3 axis,
  double3 halfAx,
  int atomCount,
  double t_max,
  ulong step,
  ulong seed,
  double BETA,
  int *kill);

#endif
