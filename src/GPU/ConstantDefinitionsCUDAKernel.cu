#include "ConstantDefinitionsCUDAKernel.cuh"

__constant__ double gpu_sigmaSq[1000];
__constant__ double gpu_epsilon_Cn[1000];
__constant__ double gpu_n[1000];
__constant__ int gpu_VDW_Kind;
__constant__ int gpu_count;
__constant__ bool gpu_isMartini;
__constant__ double gpu_rCut;
__constant__ double gpu_rCutLow;
__constant__ double gpu_rOn;
__constant__ double gpu_alpha;
__constant__ bool gpu_ewald;
__constant__ double gpu_diElectric_1;
