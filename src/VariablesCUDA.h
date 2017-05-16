#pragma once
#ifdef GOMC_CUDA
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
    gpu_rCutLow = NULL;
    gpu_rOn = NULL;
    gpu_alpha = NULL;
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
  double *gpu_rCutLow;
  double *gpu_rOn;
  double *gpu_alpha;
  int *gpu_ewald;
  double *gpu_diElectric_1;
};
#endif
