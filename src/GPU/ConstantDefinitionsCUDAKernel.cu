/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CUDAMemoryManager.cuh"
#include <iostream>
#include <stdio.h>

void InitGPULambda(VariablesCUDA *vars, int *molIndex, int *kindIndex,
                   double *lambdaVDW, double *lambdaCoulomb, bool *isFraction)
{
  // allocate gpu memory for lambda variables
  CUDAMemoryManager::mallocMemory(&vars->gpu_molIndex, (int)BOX_TOTAL * sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars->gpu_kindIndex, (int)BOX_TOTAL * sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars->gpu_lambdaVDW, (int)BOX_TOTAL * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_lambdaCoulomb, (int)BOX_TOTAL * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_isFraction, (int)BOX_TOTAL * sizeof(bool));

  // copy required data
  cudaMemcpy(vars->gpu_molIndex, molIndex, BOX_TOTAL * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_kindIndex, kindIndex, BOX_TOTAL * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_lambdaVDW, lambdaVDW, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_lambdaCoulomb, lambdaCoulomb, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_isFraction, isFraction, BOX_TOTAL * sizeof(bool),
             cudaMemcpyHostToDevice);
}

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn,
                       double const *n, int VDW_Kind, int isMartini,
                       int count, double Rcut, double const *rCutCoulomb,
                       double RcutLow, double Ron, double const *alpha,
                       int ewald, double diElectric_1)
{
  int countSq = count * count;
  CUDAMemoryManager::mallocMemory(&vars.gpu_sigmaSq, countSq * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_epsilon_Cn, countSq * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_n, countSq * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_VDW_Kind, sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars.gpu_isMartini, sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars.gpu_count, sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars.gpu_rCut, sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_rCutCoulomb, BOX_TOTAL * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_rCutLow, sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_rOn, sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_alpha, BOX_TOTAL * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars.gpu_ewald, sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars.gpu_diElectric_1, sizeof(double));

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
  cudaMemcpy(vars.gpu_rCutCoulomb, rCutCoulomb, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutLow, &RcutLow, sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rOn, &Ron, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_alpha, alpha, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_ewald, &ewald, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_diElectric_1, &diElectric_1, sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
                         uint maxAtomsInMol, uint maxMolNumber)
{
  CUDAMemoryManager::mallocMemory(&vars->gpu_x, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_y, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_z, atomNumber * sizeof(double));

  CUDAMemoryManager::mallocMemory(&vars->gpu_dx, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_dy, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_dz, atomNumber * sizeof(double));

  CUDAMemoryManager::mallocMemory(&vars->gpu_nx, maxAtomsInMol * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_ny, maxAtomsInMol * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_nz, maxAtomsInMol * sizeof(double));

  CUDAMemoryManager::mallocMemory(&vars->gpu_comx, maxMolNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_comy, maxMolNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_comz, maxMolNumber * sizeof(double));

  CUDAMemoryManager::mallocMemory(&vars->gpu_rT11, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rT12, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rT13, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rT22, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rT23, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rT33, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT11, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT12, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT13, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT22, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT23, MAX_PAIR_SIZE * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_vT33, MAX_PAIR_SIZE * sizeof(double));

  CUDAMemoryManager::mallocMemory(&vars->gpu_nonOrth, sizeof(int));
  vars->gpu_cell_x = new double *[BOX_TOTAL];
  vars->gpu_cell_y = new double *[BOX_TOTAL];
  vars->gpu_cell_z = new double *[BOX_TOTAL];
  vars->gpu_Invcell_x = new double *[BOX_TOTAL];
  vars->gpu_Invcell_y = new double *[BOX_TOTAL];
  vars->gpu_Invcell_z = new double *[BOX_TOTAL];
  for(uint b = 0; b < BOX_TOTAL; b++) {
    CUDAMemoryManager::mallocMemory(&vars->gpu_cell_x[b], 3 * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_cell_y[b], 3 * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_cell_z[b], 3 * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_Invcell_x[b], 3 * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_Invcell_y[b], 3 * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_Invcell_z[b], 3 * sizeof(double));
  }

  CUDAMemoryManager::mallocMemory(&vars->gpu_aForcex, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_aForcey, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_aForcez, atomNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_mForcex, maxMolNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_mForcey, maxMolNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_mForcez, maxMolNumber * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_cellVector, atomNumber * sizeof(int));
  CUDAMemoryManager::mallocMemory(&vars->gpu_mapParticleToCell, atomNumber * sizeof(int));
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void InitExp6Variables(VariablesCUDA *vars, double *rMin, double *expConst,
                       double *rMaxSq, uint size)
{
  CUDAMemoryManager::mallocMemory(&vars->gpu_rMin, size * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_rMaxSq, size * sizeof(double));
  CUDAMemoryManager::mallocMemory(&vars->gpu_expConst, size * sizeof(double));

  cudaMemcpy(vars->gpu_rMin, rMin, size * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_rMaxSq, rMaxSq, size * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_expConst, expConst, size * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
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

  for(uint b = 0; b < BOX_TOTAL; b++) {
    CUDAMemoryManager::mallocMemory(&vars->gpu_kx[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_ky[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_kz[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_kxRef[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_kyRef[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_kzRef[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_sumRnew[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_sumRref[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_sumInew[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_sumIref[b], imageTotal * sizeof(double));

    CUDAMemoryManager::mallocMemory(&vars->gpu_prefact[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_prefactRef[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_hsqr[b], imageTotal * sizeof(double));
    CUDAMemoryManager::mallocMemory(&vars->gpu_hsqrRef[b], imageTotal * sizeof(double));
  }
  checkLastErrorCUDA(__FILE__, __LINE__);
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
  checkLastErrorCUDA(__FILE__, __LINE__);
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

void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z)
{
  int nonOrth = 0;
  cudaMemcpy(vars->gpu_cell_x[box], cellBasis_x, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cell_y[box], cellBasis_y, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cell_z[box], cellBasis_z, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nonOrth, &nonOrth, sizeof(int), cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void UpdateInvCellBasisCUDA(VariablesCUDA *vars, uint box,
                            double *invCellBasis_x, double *invCellBasis_y,
                            double *invCellBasis_z)
{
  int nonOrth = 1;
  cudaMemcpy(vars->gpu_Invcell_x[box], invCellBasis_x, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_Invcell_y[box], invCellBasis_y, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_Invcell_z[box], invCellBasis_z, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nonOrth, &nonOrth, sizeof(int), cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void DestroyEwaldCUDAVars(VariablesCUDA *vars)
{
  for(uint b = 0; b < BOX_TOTAL; b++) {
    CUDAMemoryManager::freeMemory(vars->gpu_kx[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_ky[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_kz[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_kxRef[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_kyRef[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_kzRef[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_sumRnew[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_sumRref[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_sumInew[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_sumIref[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_prefact[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_prefactRef[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_hsqr[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_hsqrRef[b]);
  }
  delete [] vars->gpu_kx;
  delete [] vars->gpu_ky;
  delete [] vars->gpu_kz;
  delete [] vars->gpu_kxRef;
  delete [] vars->gpu_kyRef;
  delete [] vars->gpu_kzRef;
  delete [] vars->gpu_sumRnew;
  delete [] vars->gpu_sumRref;
  delete [] vars->gpu_sumInew;
  delete [] vars->gpu_sumIref;
  delete [] vars->gpu_prefact;
  delete [] vars->gpu_prefactRef;
  delete [] vars->gpu_hsqr;
  delete [] vars->gpu_hsqrRef;
}

void DestroyCUDAVars(VariablesCUDA *vars)
{
  CUDAMemoryManager::freeMemory(vars->gpu_sigmaSq);
  CUDAMemoryManager::freeMemory(vars->gpu_epsilon_Cn);
  CUDAMemoryManager::freeMemory(vars->gpu_n);
  CUDAMemoryManager::freeMemory(vars->gpu_VDW_Kind);
  CUDAMemoryManager::freeMemory(vars->gpu_isMartini);
  CUDAMemoryManager::freeMemory(vars->gpu_count);
  CUDAMemoryManager::freeMemory(vars->gpu_rCut);
  CUDAMemoryManager::freeMemory(vars->gpu_rCutCoulomb);
  CUDAMemoryManager::freeMemory(vars->gpu_rCutLow);
  CUDAMemoryManager::freeMemory(vars->gpu_rOn);
  CUDAMemoryManager::freeMemory(vars->gpu_alpha);
  CUDAMemoryManager::freeMemory(vars->gpu_ewald);
  CUDAMemoryManager::freeMemory(vars->gpu_diElectric_1);
  CUDAMemoryManager::freeMemory(vars->gpu_x);
  CUDAMemoryManager::freeMemory(vars->gpu_y);
  CUDAMemoryManager::freeMemory(vars->gpu_z);
  CUDAMemoryManager::freeMemory(vars->gpu_dx);
  CUDAMemoryManager::freeMemory(vars->gpu_dy);
  CUDAMemoryManager::freeMemory(vars->gpu_dz);
  CUDAMemoryManager::freeMemory(vars->gpu_nx);
  CUDAMemoryManager::freeMemory(vars->gpu_ny);
  CUDAMemoryManager::freeMemory(vars->gpu_nz);
  CUDAMemoryManager::freeMemory(vars->gpu_comx);
  CUDAMemoryManager::freeMemory(vars->gpu_comy);
  CUDAMemoryManager::freeMemory(vars->gpu_comz);
  CUDAMemoryManager::freeMemory(vars->gpu_aForcex);
  CUDAMemoryManager::freeMemory(vars->gpu_aForcey);
  CUDAMemoryManager::freeMemory(vars->gpu_aForcez);
  CUDAMemoryManager::freeMemory(vars->gpu_mForcex);
  CUDAMemoryManager::freeMemory(vars->gpu_mForcey);
  CUDAMemoryManager::freeMemory(vars->gpu_mForcez);
  CUDAMemoryManager::freeMemory(vars->gpu_cellVector);
  CUDAMemoryManager::freeMemory(vars->gpu_mapParticleToCell);
  CUDAMemoryManager::freeMemory(vars->gpu_rT11);
  CUDAMemoryManager::freeMemory(vars->gpu_rT12);
  CUDAMemoryManager::freeMemory(vars->gpu_rT13);
  CUDAMemoryManager::freeMemory(vars->gpu_rT22);
  CUDAMemoryManager::freeMemory(vars->gpu_rT23);
  CUDAMemoryManager::freeMemory(vars->gpu_rT33);
  CUDAMemoryManager::freeMemory(vars->gpu_vT11);
  CUDAMemoryManager::freeMemory(vars->gpu_vT12);
  CUDAMemoryManager::freeMemory(vars->gpu_vT13);
  CUDAMemoryManager::freeMemory(vars->gpu_vT22);
  CUDAMemoryManager::freeMemory(vars->gpu_vT23);
  CUDAMemoryManager::freeMemory(vars->gpu_vT33);
  CUDAMemoryManager::freeMemory(vars->gpu_nonOrth);
  for(uint b = 0; b < BOX_TOTAL; b++) {
    CUDAMemoryManager::freeMemory(vars->gpu_cell_x[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_cell_y[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_cell_z[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_Invcell_x[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_Invcell_y[b]);
    CUDAMemoryManager::freeMemory(vars->gpu_Invcell_z[b]);
  }
  delete [] vars-> gpu_cell_x;
  delete [] vars-> gpu_cell_y;
  delete [] vars-> gpu_cell_z;
  delete [] vars-> gpu_Invcell_x;
  delete [] vars-> gpu_Invcell_y;
  delete [] vars-> gpu_Invcell_z;
}

#endif /*GOMC_CUDA*/
