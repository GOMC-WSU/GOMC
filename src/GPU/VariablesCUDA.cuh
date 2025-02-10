/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include "EnsemblePreprocessor.h"
#include "NumLib.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

// Need a separate float constant for device code with the MSVC compiler
// See CUDA Programming Guide section I.4.13 for details
static const __device__ double qqFactGPU = num::qqFact;

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

inline void checkLastErrorCUDA(const char *file, int line) {
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    exit(code);
  }
}

inline void printFreeMemory() {
  size_t free_byte;
  size_t total_byte;
  cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

  if (cudaSuccess != cuda_status) {
    printf("Error: cudaMemGetInfo fails, %s \n",
           cudaGetErrorString(cuda_status));
    exit(1);
  }
  double free_db = (double)free_byte;
  double total_db = (double)total_byte;
  double used_db = total_db - free_db;
  printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
         used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0,
         total_db / 1024.0 / 1024.0);
}

class VariablesCUDA {
public:
  VariablesCUDA() {
    gpu_sigmaSq = NULL;
    gpu_epsilon_Cn = NULL;
    gpu_n = NULL;
    gpu_VDW_Kind = NULL;
    gpu_isMartini = NULL;
    gpu_count = NULL;
    gpu_rCut = NULL;
    gpu_rCutSq = NULL;
    gpu_rCutLow = NULL;
    gpu_rOn = NULL;
    gpu_alpha = NULL;
    gpu_alphaSq = NULL;
    gpu_rCutCoulomb = NULL;
    gpu_rCutCoulombSq = NULL;
    gpu_ewald = NULL;
    gpu_diElectric_1 = NULL;
    gpu_aForcex = NULL;
    gpu_aForcey = NULL;
    gpu_aForcez = NULL;
    gpu_mForcex = NULL;
    gpu_mForcey = NULL;
    gpu_mForcez = NULL;
    gpu_startAtomIdx = NULL;

    // setting lambda values to null
    gpu_molIndex = NULL;
    gpu_lambdaVDW = NULL;
    gpu_lambdaCoulomb = NULL;
    gpu_isFraction = NULL;
  }
  double *gpu_sigmaSq;
  double *gpu_epsilon_Cn;
  double *gpu_n;
  int *gpu_VDW_Kind;
  int *gpu_isMartini;
  int *gpu_count;
  int *gpu_startAtomIdx; // start atom index of the molecule
  double *gpu_rCut, *gpu_rCutSq;
  double *gpu_rCutCoulomb, *gpu_rCutCoulombSq;
  double *gpu_rCutLow;
  double *gpu_rOn;
  double *gpu_alpha, *gpu_alphaSq;
  int *gpu_ewald;
  double *gpu_diElectric_1;
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_nx, *gpu_ny, *gpu_nz;
  double *gpu_dx, *gpu_dy, *gpu_dz;
  double **gpu_kx, **gpu_ky, **gpu_kz;
  double **gpu_kxRef, **gpu_kyRef, **gpu_kzRef;
  double **gpu_sumRnew, **gpu_sumInew, **gpu_sumRref, **gpu_sumIref;
  double **gpu_prefact, **gpu_prefactRef;
  double **gpu_hsqr, **gpu_hsqrRef;
  double *gpu_comx, *gpu_comy, *gpu_comz;
  double *gpu_rT11, *gpu_rT12, *gpu_rT13;
  double *gpu_rT22, *gpu_rT23, *gpu_rT33;
  double *gpu_vT11, *gpu_vT12, *gpu_vT13;
  double *gpu_vT22, *gpu_vT23, *gpu_vT33;
  double **gpu_cell_x, **gpu_cell_y, **gpu_cell_z;
  double **gpu_Invcell_x, **gpu_Invcell_y, **gpu_Invcell_z;
  int *gpu_nonOrth;
  double *gpu_aForcex, *gpu_aForcey, *gpu_aForcez;
  double *gpu_mForcex, *gpu_mForcey, *gpu_mForcez;
  double *gpu_mTorquex, *gpu_mTorquey, *gpu_mTorquez;
  int *gpu_inForceRange;
  double *gpu_aForceRecx, *gpu_aForceRecy, *gpu_aForceRecz;
  double *gpu_mForceRecx, *gpu_mForceRecy, *gpu_mForceRecz;
  double *gpu_rMin, *gpu_expConst, *gpu_rMaxSq;

  double *gpu_r_k_x, *gpu_r_k_y, *gpu_r_k_z;
  double *gpu_t_k_x, *gpu_t_k_y, *gpu_t_k_z;

  // lambda structure
  int *gpu_molIndex;
  double *gpu_lambdaVDW, *gpu_lambdaCoulomb;
  bool *gpu_isFraction;

  // new pair interaction calculation done on GPU
  int *gpu_cellVector, *gpu_mapParticleToCell;
};
#endif
