/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include "EnsemblePreprocessor.h"
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
  }
  double *gpu_sigmaSq;
  double *gpu_epsilon_Cn;
  double *gpu_n;
  int *gpu_VDW_Kind;
  int *gpu_isMartini;
  int *gpu_count;
  double *gpu_rCut;
  double *gpu_rCutCoulomb;
  double *gpu_rCutLow;
  double *gpu_rOn;
  double *gpu_alpha;
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
};
#endif
