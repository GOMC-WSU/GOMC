/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include <cuda.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "EnsemblePreprocessor.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

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
  real free_db = (real)free_byte ;
  real total_db = (real)total_byte ;
  real used_db = total_db - free_db ;
  printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
         used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
}

class VariablesCUDA
{
public:
  VariablesCUDA()
  {
    gpu_sigmaSq = NULL;
    gpu_epsilon_Cn = NULL;
    gpu_n = NULL;
    gpu_VDW_Kind = NULL;
    gpu_isMartini = NULL;
    gpu_count = NULL;
    gpu_rCut = NULL;
    gpu_rCutLow = NULL;
    gpu_rOn = NULL;
    gpu_alpha = NULL;
    gpu_rCutCoulomb = NULL;
    gpu_ewald = NULL;
    gpu_diElectric_1 = NULL;
    gpu_aForcex = NULL;
    gpu_aForcey = NULL;
    gpu_aForcez = NULL;
    gpu_mForcex = NULL;
    gpu_mForcey = NULL;
    gpu_mForcez = NULL;
  }
  real *gpu_sigmaSq;
  real *gpu_epsilon_Cn;
  real *gpu_n;
  int *gpu_VDW_Kind;
  int *gpu_isMartini;
  int *gpu_count;
  real *gpu_rCut;
  real *gpu_rCutCoulomb;
  real *gpu_rCutLow;
  real *gpu_rOn;
  real *gpu_alpha;
  int *gpu_ewald;
  real *gpu_diElectric_1;
  real *gpu_x, *gpu_y, *gpu_z;
  real *gpu_nx, *gpu_ny, *gpu_nz;
  real *gpu_dx, *gpu_dy, *gpu_dz;
  real **gpu_kx, **gpu_ky, **gpu_kz;
  real **gpu_kxRef, **gpu_kyRef, **gpu_kzRef;
  real **gpu_sumRnew, **gpu_sumInew, **gpu_sumRref, **gpu_sumIref;
  real **gpu_prefact, **gpu_prefactRef;
  real **gpu_hsqr, **gpu_hsqrRef;
  real *gpu_comx, *gpu_comy, *gpu_comz;
  real *gpu_rT11, *gpu_rT12, *gpu_rT13;
  real *gpu_rT22, *gpu_rT23, *gpu_rT33;
  real *gpu_vT11, *gpu_vT12, *gpu_vT13;
  real *gpu_vT22, *gpu_vT23, *gpu_vT33;
  real **gpu_cell_x, **gpu_cell_y, **gpu_cell_z;
  real **gpu_Invcell_x, **gpu_Invcell_y, **gpu_Invcell_z;
  int *gpu_nonOrth;
  real *gpu_aForcex, *gpu_aForcey, *gpu_aForcez;
  real *gpu_mForcex, *gpu_mForcey, *gpu_mForcez;
};
#endif
