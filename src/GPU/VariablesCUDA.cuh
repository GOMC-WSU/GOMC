/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "EnsemblePreprocessor.h"
#include "NumLib.h"

//Need a separate float constant for device code with the MSVC compiler
//See CUDA Programming Guide section I.4.13 for details 
static const __device__ double qqFactGPU = num::qqFact;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

#ifndef NDEBUG
inline void checkLastErrorCUDA(const char *file, int line)
{
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    exit(code);
  }
}
#endif

inline void printFreeMemory()
{
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

  if ( cudaSuccess != cuda_status ) {
    printf("Error: cudaMemGetInfo fails, %s \n",
           cudaGetErrorString(cuda_status) );
    exit(1);
  }
  double free_db = (double)free_byte ;
  double total_db = (double)total_byte ;
  double used_db = total_db - free_db ;
  printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
         used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
}

class VariablesCUDA
{
public:
  VariablesCUDA()
  {
    cub_reduce_storage_size = 0;
    cub_reduce_storage = nullptr;
    gpu_sigmaSq = nullptr;
    gpu_epsilon_Cn = nullptr;
    gpu_n = nullptr;
    gpu_VDW_Kind = nullptr;
    gpu_isMartini = nullptr;
    gpu_count = nullptr;
    gpu_startAtomIdx = nullptr;
    gpu_rCut = nullptr;
    gpu_rCutCoulomb = nullptr;
    gpu_rCutLow = nullptr;
    gpu_rOn = nullptr;
    gpu_alpha = nullptr;
    gpu_ewald = nullptr;
    gpu_diElectric_1 = nullptr;
    gpu_finalVal = nullptr;
    gpu_aForcex = nullptr;
    gpu_aForcey = nullptr;
    gpu_aForcez = nullptr;
    gpu_mForcex = nullptr;
    gpu_mForcey = nullptr;
    gpu_mForcez = nullptr;
    gpu_startAtomIdx = nullptr;

    // setting lambda valuesy to null
    gpu_molIndex = nullptr;
    gpu_lambdaVDW = nullptr;
    gpu_lambdaCoulomb = nullptr;
    gpu_isFraction = nullptr;
  }

  size_t cub_reduce_storage_size;
  void *cub_reduce_storage;
  double *gpu_sigmaSq;
  double *gpu_epsilon_Cn;
  double *gpu_n;
  int *gpu_VDW_Kind;
  int *gpu_isMartini;
  int *gpu_count;
  int *gpu_startAtomIdx; //start atom index of the molecule
  double *gpu_rCut;
  double *gpu_rCutCoulomb;
  double *gpu_rCutLow;
  double *gpu_rOn;
  double *gpu_alpha;
  int *gpu_ewald;
  double *gpu_diElectric_1;
  double *gpu_finalVal;
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_nx, *gpu_ny, *gpu_nz;
  double *gpu_dx, *gpu_dy, *gpu_dz;
  double **gpu_kx, **gpu_ky, **gpu_kz;
  double **gpu_kxRef, **gpu_kyRef, **gpu_kzRef;
  double **gpu_sumRnew, **gpu_sumInew, **gpu_sumRref, **gpu_sumIref;
  double **gpu_prefact, **gpu_prefactRef;
  double **gpu_hsqr, **gpu_hsqrRef;
  double *gpu_recipEnergies;
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
