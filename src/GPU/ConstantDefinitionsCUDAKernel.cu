#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include <iostream>

#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq, 
		       double const *epsilon_Cn,
		       double const *n, int VDW_Kind, int isMartini,
		       int count, double Rcut, double RcutLow,
		       double Ron, double alpha,
		       int ewald, double diElectric_1)
{
  int countSq = count * count;
  cudaMalloc(&vars.gpu_sigmaSq, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_epsilon_Cn, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_n, countSq * sizeof(double));
  cudaMalloc(&vars.gpu_VDW_Kind, sizeof(int));
  cudaMalloc(&vars.gpu_isMartini, sizeof(int));
  cudaMalloc(&vars.gpu_count, sizeof(int));
  cudaMalloc(&vars.gpu_rCut, sizeof(double));
  cudaMalloc(&vars.gpu_rCutLow, sizeof(double));
  cudaMalloc(&vars.gpu_rOn, sizeof(double));
  cudaMalloc(&vars.gpu_alpha, sizeof(double));
  cudaMalloc(&vars.gpu_ewald, sizeof(int));
  cudaMalloc(&vars.gpu_diElectric_1, sizeof(double));

  cudaMemcpy(vars.gpu_sigmaSq, sigmaSq, countSq * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_epsilon_Cn, epsilon_Cn, countSq * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_n, n, countSq * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_VDW_Kind, &VDW_Kind, sizeof(int),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_isMartini, &isMartini, sizeof(int),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_count, &count, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCut, &Rcut, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutLow, &RcutLow, sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rOn, &Ron, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_alpha, &alpha, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_ewald, &ewald, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_diElectric_1, &diElectric_1, sizeof(double),
	     cudaMemcpyHostToDevice);
}

void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
			 uint maxAtomsInMol, uint maxMolNumber)
{
  uint maxPair = (atomNumber * atomNumber) / 2;
  cudaMalloc(&vars->gpu_x, atomNumber * sizeof(double));
  cudaMalloc(&vars->gpu_y, atomNumber * sizeof(double));
  cudaMalloc(&vars->gpu_z, atomNumber * sizeof(double));
  
  cudaMalloc(&vars->gpu_nx, maxAtomsInMol * sizeof(double));
  cudaMalloc(&vars->gpu_ny, maxAtomsInMol * sizeof(double));
  cudaMalloc(&vars->gpu_nz, maxAtomsInMol * sizeof(double));
  
  cudaMalloc(&vars->gpu_comx, maxMolNumber * sizeof(double));
  cudaMalloc(&vars->gpu_comy, maxMolNumber * sizeof(double));
  cudaMalloc(&vars->gpu_comz, maxMolNumber * sizeof(double));

  cudaMalloc(&vars->gpu_rT11, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_rT12, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_rT13, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_rT22, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_rT23, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_rT33, maxPair * sizeof(double));

  cudaMalloc(&vars->gpu_vT11, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_vT12, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_vT13, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_vT22, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_vT23, maxPair * sizeof(double));
  cudaMalloc(&vars->gpu_vT33, maxPair * sizeof(double));
}

void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal)
{
  vars->gpu_kx = new double *[BOX_TOTAL];
  vars->gpu_ky = new double *[BOX_TOTAL];
  vars->gpu_kz = new double *[BOX_TOTAL];
  vars->gpu_kxRef = new double *[BOX_TOTAL];
  vars->gpu_kyRef = new double *[BOX_TOTAL];
  vars->gpu_kzRef = new double *[BOX_TOTAL];
  vars->gpu_sumRnew = new double *[BOX_TOTAL];
  vars->gpu_sumRref = new double *[BOX_TOTAL];
  vars->gpu_sumInew = new double *[BOX_TOTAL];
  vars->gpu_sumIref = new double *[BOX_TOTAL];
  vars->gpu_prefact = new double *[BOX_TOTAL];
  vars->gpu_prefactRef = new double *[BOX_TOTAL];
  vars->gpu_hsqr = new double *[BOX_TOTAL];
  vars->gpu_hsqrRef = new double *[BOX_TOTAL];
  for(uint b = 0; b < BOX_TOTAL; b++)
  {
    cudaMalloc(&vars->gpu_kx[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_ky[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_kz[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_kxRef[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_kyRef[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_kzRef[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_sumRnew[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_sumRref[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_sumInew[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_sumIref[b], imageTotal * sizeof(double));
    
    cudaMalloc(&vars->gpu_prefact[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_prefactRef[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_hsqr[b], imageTotal * sizeof(double));
    cudaMalloc(&vars->gpu_hsqrRef[b], imageTotal * sizeof(double));
  }
}

void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal)
{
  cudaMemcpy(vars->gpu_sumRref[box], vars->gpu_sumRnew[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_sumIref[box], vars->gpu_sumInew[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_prefactRef[box], vars->gpu_prefact[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_hsqrRef[box], vars->gpu_hsqr[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_kxRef[box], vars->gpu_kx[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_kyRef[box], vars->gpu_ky[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_kzRef[box], vars->gpu_kz[box],
	     imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
}

void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box)
{
  double *tempKx, *tempKy, *tempKz, *tempHsqr, *tempPrefact;
  tempKx = vars->gpu_kxRef[box];
  tempKy = vars->gpu_kyRef[box];
  tempKz = vars->gpu_kzRef[box];
  tempHsqr = vars->gpu_hsqrRef[box];
  tempPrefact = vars->gpu_prefactRef[box];
  
  vars->gpu_kxRef[box] = vars->gpu_kx[box];
  vars->gpu_kyRef[box] = vars->gpu_ky[box];
  vars->gpu_kzRef[box] = vars->gpu_kz[box];
  vars->gpu_hsqrRef[box] = vars->gpu_hsqr[box];
  vars->gpu_prefactRef[box] = vars->gpu_prefact[box];
  
  vars->gpu_kx[box] = tempKx;
  vars->gpu_ky[box] = tempKy;
  vars->gpu_kz[box] = tempKz;
  vars->gpu_hsqr[box] = tempHsqr;
  vars->gpu_prefact[box] = tempPrefact;
}

void UpdateRecipCUDA(VariablesCUDA *vars, uint box)
{
  double *tempR, *tempI;
  tempR = vars->gpu_sumRref[box];
  tempI = vars->gpu_sumIref[box];
  vars->gpu_sumRref[box] = vars->gpu_sumRnew[box];
  vars->gpu_sumIref[box] = vars->gpu_sumInew[box];
  vars->gpu_sumRnew[box] = tempR;
  vars->gpu_sumInew[box] = tempI;
}

#endif /*GOMC_CUDA*/
