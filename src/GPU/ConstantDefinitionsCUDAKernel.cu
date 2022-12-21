/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include <iostream>

#include "CUDAMemoryManager.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "GeomLib.h"

void UpdateGPULambda(VariablesCUDA *vars, int *molIndex, double *lambdaVDW,
                     double *lambdaCoulomb, bool *isFraction) {
  // copy lambda data
  cudaMemcpy(vars->gpu_molIndex, molIndex, BOX_TOTAL * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_lambdaVDW, lambdaVDW, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_lambdaCoulomb, lambdaCoulomb, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_isFraction, isFraction, BOX_TOTAL * sizeof(bool),
             cudaMemcpyHostToDevice);
}

void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n, int VDW_Kind,
                       int isMartini, int count, double Rcut,
                       double const *rCutCoulomb, double RcutLow, double Ron,
                       double const *alpha, int ewald, double diElectric_1,
                       int wolf, 
                       int coulKind,
                       double const * wolfAlpha,
                       double const * wolfFactor1, 
                       double const * wolfFactor2, 
                       double const * wolfFactor3) {
  int countSq = count * count;
  CUMALLOC((void **)&vars.gpu_sigmaSq, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_epsilon_Cn, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_n, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_VDW_Kind, sizeof(int));
  CUMALLOC((void **)&vars.gpu_isMartini, sizeof(int));
  CUMALLOC((void **)&vars.gpu_count, sizeof(int));
  CUMALLOC((void **)&vars.gpu_rCut, sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutCoulomb, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutLow, sizeof(double));
  CUMALLOC((void **)&vars.gpu_rOn, sizeof(double));
  CUMALLOC((void **)&vars.gpu_alpha, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_ewald, sizeof(int));
  CUMALLOC((void **)&vars.gpu_diElectric_1, sizeof(double));

  // allocate gpu memory for lambda variables
  CUMALLOC((void **)&vars.gpu_molIndex, (int)BOX_TOTAL * sizeof(int));
  CUMALLOC((void **)&vars.gpu_lambdaVDW, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_lambdaCoulomb, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_isFraction, (int)BOX_TOTAL * sizeof(bool));


  // allocate gpu memory for wolf variables
  CUMALLOC((void**) &vars.gpu_wolf, sizeof(int));
  CUMALLOC((void**) &vars.gpu_coulKind, sizeof(int));
  CUMALLOC((void**) &vars.gpu_wolfAlpha, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void**) &vars.gpu_wolfFactor1, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void**) &vars.gpu_wolfFactor2, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void**) &vars.gpu_wolfFactor3, (int)BOX_TOTAL * sizeof(double));

  cudaMemcpy(vars.gpu_sigmaSq, sigmaSq, countSq * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_epsilon_Cn, epsilon_Cn, countSq * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_n, n, countSq * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_VDW_Kind, &VDW_Kind, sizeof(int), cudaMemcpyHostToDevice);
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

  cudaMemcpy(vars.gpu_wolf, &wolf, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_coulKind, &coulKind, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfAlpha, wolfAlpha, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfFactor1, wolfFactor1, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfFactor2, wolfFactor2, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice); 
  cudaMemcpy(vars.gpu_wolfFactor3, wolfFactor3, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);  
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
                         uint maxAtomsInMol, uint maxMolNumber) {
  CUMALLOC((void **)&vars->gpu_x, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_y, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_z, atomNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_dx, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_dy, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_dz, atomNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_nx, maxAtomsInMol * sizeof(double));
  CUMALLOC((void **)&vars->gpu_ny, maxAtomsInMol * sizeof(double));
  CUMALLOC((void **)&vars->gpu_nz, maxAtomsInMol * sizeof(double));

  CUMALLOC((void **)&vars->gpu_comx, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_comy, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_comz, maxMolNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_r_k_x, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_r_k_y, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_r_k_z, maxMolNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_t_k_x, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_t_k_y, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_t_k_z, maxMolNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_nonOrth, sizeof(int));
  vars->gpu_cell_x = new double *[BOX_TOTAL];
  vars->gpu_cell_y = new double *[BOX_TOTAL];
  vars->gpu_cell_z = new double *[BOX_TOTAL];
  vars->gpu_Invcell_x = new double *[BOX_TOTAL];
  vars->gpu_Invcell_y = new double *[BOX_TOTAL];
  vars->gpu_Invcell_z = new double *[BOX_TOTAL];
  for (uint b = 0; b < BOX_TOTAL; b++) {
    CUMALLOC((void **)&vars->gpu_cell_x[b], 3 * sizeof(double));
    CUMALLOC((void **)&vars->gpu_cell_y[b], 3 * sizeof(double));
    CUMALLOC((void **)&vars->gpu_cell_z[b], 3 * sizeof(double));
    CUMALLOC((void **)&vars->gpu_Invcell_x[b], 3 * sizeof(double));
    CUMALLOC((void **)&vars->gpu_Invcell_y[b], 3 * sizeof(double));
    CUMALLOC((void **)&vars->gpu_Invcell_z[b], 3 * sizeof(double));
  }

  CUMALLOC((void **)&vars->gpu_aForcex, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForcey, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForcez, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcex, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcey, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcez, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquex, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquey, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquez, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_inForceRange, maxMolNumber * sizeof(int));
  CUMALLOC((void **)&vars->gpu_aForceRecx, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForceRecy, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForceRecz, atomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecx, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecy, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecz, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_cellVector, atomNumber * sizeof(int));
  CUMALLOC((void **)&vars->gpu_mapParticleToCell, atomNumber * sizeof(int));
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void UpdateGPUWolfEwald(VariablesCUDA &vars, 
                       int ewald,
                       int wolf, 
                       int coulKind,
                       double const * wolfAlpha,
                       double const * wolfFactor1, 
                       double const * wolfFactor2, 
                       double const * wolfFactor3){

  cudaMemcpy(vars.gpu_ewald, &ewald, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolf, &wolf, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_coulKind, &coulKind, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfAlpha, wolfAlpha, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfFactor1, wolfFactor1, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_wolfFactor2, wolfFactor2, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice); 
  cudaMemcpy(vars.gpu_wolfFactor3, wolfFactor3, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);  
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void InitExp6Variables(VariablesCUDA *vars, double *rMin, double *expConst,
                       double *rMaxSq, uint size) {
  CUMALLOC((void **)&vars->gpu_rMin, size * sizeof(double));
  CUMALLOC((void **)&vars->gpu_rMaxSq, size * sizeof(double));
  CUMALLOC((void **)&vars->gpu_expConst, size * sizeof(double));

  cudaMemcpy(vars->gpu_rMin, rMin, size * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_rMaxSq, rMaxSq, size * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_expConst, expConst, size * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal) {
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

  for (uint b = 0; b < BOX_TOTAL; b++) {
    CUMALLOC((void **)&vars->gpu_kx[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_ky[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_kz[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_kxRef[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_kyRef[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_kzRef[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_sumRnew[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_sumRref[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_sumInew[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_sumIref[b], imageTotal * sizeof(double));

    CUMALLOC((void **)&vars->gpu_prefact[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_prefactRef[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_hsqr[b], imageTotal * sizeof(double));
    CUMALLOC((void **)&vars->gpu_hsqrRef[b], imageTotal * sizeof(double));
  }
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal) {
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

void CopyRefToNewCUDA(VariablesCUDA *vars, uint box, uint imageTotal) {
  cudaMemcpy(vars->gpu_sumRnew[box], vars->gpu_sumRref[box],
             imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_sumInew[box], vars->gpu_sumIref[box],
             imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box) {
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

void UpdateRecipCUDA(VariablesCUDA *vars, uint box) {
  double *tempR, *tempI;
  tempR = vars->gpu_sumRref[box];
  tempI = vars->gpu_sumIref[box];
  vars->gpu_sumRref[box] = vars->gpu_sumRnew[box];
  vars->gpu_sumIref[box] = vars->gpu_sumInew[box];
  vars->gpu_sumRnew[box] = tempR;
  vars->gpu_sumInew[box] = tempI;
}

void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z) {
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
                            double *invCellBasis_z) {
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

void DestroyEwaldCUDAVars(VariablesCUDA *vars) {
  for (uint b = 0; b < BOX_TOTAL; b++) {
    CUFREE(vars->gpu_kx[b]);
    CUFREE(vars->gpu_ky[b]);
    CUFREE(vars->gpu_kz[b]);
    CUFREE(vars->gpu_kxRef[b]);
    CUFREE(vars->gpu_kyRef[b]);
    CUFREE(vars->gpu_kzRef[b]);
    CUFREE(vars->gpu_sumRnew[b]);
    CUFREE(vars->gpu_sumRref[b]);
    CUFREE(vars->gpu_sumInew[b]);
    CUFREE(vars->gpu_sumIref[b]);
    CUFREE(vars->gpu_prefact[b]);
    CUFREE(vars->gpu_prefactRef[b]);
    CUFREE(vars->gpu_hsqr[b]);
    CUFREE(vars->gpu_hsqrRef[b]);
  }
  delete[] vars->gpu_kx;
  delete[] vars->gpu_ky;
  delete[] vars->gpu_kz;
  delete[] vars->gpu_kxRef;
  delete[] vars->gpu_kyRef;
  delete[] vars->gpu_kzRef;
  delete[] vars->gpu_sumRnew;
  delete[] vars->gpu_sumRref;
  delete[] vars->gpu_sumInew;
  delete[] vars->gpu_sumIref;
  delete[] vars->gpu_prefact;
  delete[] vars->gpu_prefactRef;
  delete[] vars->gpu_hsqr;
  delete[] vars->gpu_hsqrRef;
}

void DestroyCUDAVars(VariablesCUDA *vars) {
  CUFREE(vars->gpu_sigmaSq);
  CUFREE(vars->gpu_epsilon_Cn);
  CUFREE(vars->gpu_n);
  CUFREE(vars->gpu_VDW_Kind);
  CUFREE(vars->gpu_isMartini);
  CUFREE(vars->gpu_count);
  CUFREE(vars->gpu_rCut);
  CUFREE(vars->gpu_rCutCoulomb);
  CUFREE(vars->gpu_rCutLow);
  CUFREE(vars->gpu_rOn);
  CUFREE(vars->gpu_alpha);
  CUFREE(vars->gpu_ewald);
  CUFREE(vars->gpu_wolf);
  CUFREE(vars->gpu_coulKind);
  CUFREE(vars->gpu_wolfAlpha);
  CUFREE(vars->gpu_wolfFactor1);
  CUFREE(vars->gpu_wolfFactor2);
  CUFREE(vars->gpu_wolfFactor3);
  CUFREE(vars->gpu_diElectric_1);
  CUFREE(vars->gpu_x);
  CUFREE(vars->gpu_y);
  CUFREE(vars->gpu_z);
  CUFREE(vars->gpu_dx);
  CUFREE(vars->gpu_dy);
  CUFREE(vars->gpu_dz);
  CUFREE(vars->gpu_nx);
  CUFREE(vars->gpu_ny);
  CUFREE(vars->gpu_nz);
  CUFREE(vars->gpu_comx);
  CUFREE(vars->gpu_comy);
  CUFREE(vars->gpu_comz);
  CUFREE(vars->gpu_r_k_x);
  CUFREE(vars->gpu_r_k_y);
  CUFREE(vars->gpu_r_k_z);
  CUFREE(vars->gpu_t_k_x);
  CUFREE(vars->gpu_t_k_y);
  CUFREE(vars->gpu_t_k_z);
  CUFREE(vars->gpu_aForcex);
  CUFREE(vars->gpu_aForcey);
  CUFREE(vars->gpu_aForcez);
  CUFREE(vars->gpu_mForcex);
  CUFREE(vars->gpu_mForcey);
  CUFREE(vars->gpu_mForcez);
  CUFREE(vars->gpu_mTorquex);
  CUFREE(vars->gpu_mTorquey);
  CUFREE(vars->gpu_mTorquez);
  CUFREE(vars->gpu_inForceRange);
  CUFREE(vars->gpu_aForceRecx);
  CUFREE(vars->gpu_aForceRecy);
  CUFREE(vars->gpu_aForceRecz);
  CUFREE(vars->gpu_mForceRecx);
  CUFREE(vars->gpu_mForceRecy);
  CUFREE(vars->gpu_mForceRecz);
  CUFREE(vars->gpu_cellVector);
  CUFREE(vars->gpu_mapParticleToCell);
  CUFREE(vars->gpu_nonOrth);
  CUFREE(vars->gpu_startAtomIdx);
  for (uint b = 0; b < BOX_TOTAL; b++) {
    CUFREE(vars->gpu_cell_x[b]);
    CUFREE(vars->gpu_cell_y[b]);
    CUFREE(vars->gpu_cell_z[b]);
    CUFREE(vars->gpu_Invcell_x[b]);
    CUFREE(vars->gpu_Invcell_y[b]);
    CUFREE(vars->gpu_Invcell_z[b]);
  }

  // delete gpu memory for lambda variables
  CUFREE(vars->gpu_molIndex);
  CUFREE(vars->gpu_lambdaVDW);
  CUFREE(vars->gpu_lambdaCoulomb);
  CUFREE(vars->gpu_isFraction);

  delete[] vars->gpu_cell_x;
  delete[] vars->gpu_cell_y;
  delete[] vars->gpu_cell_z;
  delete[] vars->gpu_Invcell_x;
  delete[] vars->gpu_Invcell_y;
  delete[] vars->gpu_Invcell_z;
}

#endif /*GOMC_CUDA*/
