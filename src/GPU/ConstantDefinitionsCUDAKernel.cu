/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA
#include "cub/cub.cuh"
#include <cuda.h>
#include <cuda_runtime.h>

#include <cstdio>
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
                       int isMartini, int count, double Rcut, double RcutSq,
                       double const *rCutCoulomb, double const *rCutCoulombSq,
                       double RcutLow, double Ron, double const *alpha,
                       double const *alphaSq, int ewald, double diElectric_1) {
  int countSq = count * count;
  CUMALLOC((void **)&vars.gpu_sigmaSq, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_epsilon_Cn, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_n, countSq * sizeof(double));
  CUMALLOC((void **)&vars.gpu_VDW_Kind, sizeof(int));
  CUMALLOC((void **)&vars.gpu_isMartini, sizeof(int));
  CUMALLOC((void **)&vars.gpu_count, sizeof(int));
  CUMALLOC((void **)&vars.gpu_rCut, sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutSq, sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutCoulomb, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutCoulombSq, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_rCutLow, sizeof(double));
  CUMALLOC((void **)&vars.gpu_rOn, sizeof(double));
  CUMALLOC((void **)&vars.gpu_alpha, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_alphaSq, BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_ewald, sizeof(int));
  CUMALLOC((void **)&vars.gpu_diElectric_1, sizeof(double));
  CUMALLOC((void **)&vars.gpu_finalVal, sizeof(double));

  // allocate GPU memory for lambda variables
  CUMALLOC((void **)&vars.gpu_molIndex, (int)BOX_TOTAL * sizeof(int));
  CUMALLOC((void **)&vars.gpu_lambdaVDW, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_lambdaCoulomb, (int)BOX_TOTAL * sizeof(double));
  CUMALLOC((void **)&vars.gpu_isFraction, (int)BOX_TOTAL * sizeof(bool));

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
  cudaMemcpy(vars.gpu_rCutSq, &RcutSq, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutCoulomb, rCutCoulomb, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutCoulombSq, rCutCoulombSq, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rCutLow, &RcutLow, sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_rOn, &Ron, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_alpha, alpha, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_alphaSq, alphaSq, BOX_TOTAL * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_ewald, &ewald, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars.gpu_diElectric_1, &diElectric_1, sizeof(double),
             cudaMemcpyHostToDevice);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

// maxAtomNumber is the maximum number of atoms in the system
// maxAtomsInMol is the maximum number of atoms in one molecule
// maxMolNumber is the maximum number of molecules in the system
void InitCoordinatesCUDA(VariablesCUDA *vars, uint maxAtomNumber,
                         uint maxAtomsInMol, uint maxMolNumber) {
  CUMALLOC((void **)&vars->gpu_x, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_y, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_z, maxAtomNumber * sizeof(double));

  // gpu_molCharge is used for a subset of the particles;
  // gpu_particleCharge is used for all the particles and is loaded only once
  CUMALLOC((void **)&vars->gpu_molCharge, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_particleCharge, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_particleKind, maxAtomNumber * sizeof(int));
  CUMALLOC((void **)&vars->gpu_particleMol, maxAtomNumber * sizeof(int));

  CUMALLOC((void **)&vars->gpu_dx, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_dy, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_dz, maxAtomNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_nx, maxAtomsInMol * sizeof(double));
  CUMALLOC((void **)&vars->gpu_ny, maxAtomsInMol * sizeof(double));
  CUMALLOC((void **)&vars->gpu_nz, maxAtomsInMol * sizeof(double));

  CUMALLOC((void **)&vars->gpu_comx, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_comy, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_comz, maxMolNumber * sizeof(double));

  CUMALLOC((void **)&vars->gpu_rt_k_x, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_rt_k_y, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_rt_k_z, maxMolNumber * sizeof(double));

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

  CUMALLOC((void **)&vars->gpu_aForcex, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForcey, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForcez, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcex, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcey, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForcez, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquex, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquey, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mTorquez, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_inForceRange, maxMolNumber * sizeof(int));
  CUMALLOC((void **)&vars->gpu_aForceRecx, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForceRecy, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_aForceRecz, maxAtomNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecx, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecy, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_mForceRecz, maxMolNumber * sizeof(double));
  CUMALLOC((void **)&vars->gpu_cellVector, maxAtomNumber * sizeof(int));
  CUMALLOC((void **)&vars->gpu_mapParticleToCell, maxAtomNumber * sizeof(int));
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

void InitExp6VariablesCUDA(VariablesCUDA *vars, double *rMin, double *expConst,
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

void InitMoleculeVariablesCUDA(VariablesCUDA *vars, const Molecules &mols) {
  const int numMol = mols.count + 1;
  // allocate memory to store molecule start atom index
  CUMALLOC((void **)&vars->gpu_startAtomIdx, numMol * sizeof(int));
  // copy start atom index
  cudaMemcpy(vars->gpu_startAtomIdx, mols.start, numMol * sizeof(int),
             cudaMemcpyHostToDevice);
}

void InitPartVariablesCUDA(VariablesCUDA *vars,
                           const std::vector<int> &particleKind,
                           const std::vector<int> &particleMol,
                           const std::vector<double> &particleCharge) {
  cudaMemcpy(vars->gpu_particleKind, &particleKind[0],
             particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_particleMol, &particleMol[0],
             particleMol.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
}

void InitEwaldVariablesCUDA(VariablesCUDA *vars, int numAtoms,
                            uint imageTotal) {
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
  CUMALLOC((void **)&vars->gpu_wT11, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_wT12, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_wT13, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_wT22, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_wT23, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_wT33, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_recipEnergies, imageTotal * sizeof(double));
  CUMALLOC((void **)&vars->gpu_particleUsed, numAtoms * sizeof(int));
  // Allocate space for cub reduction operations on the Ewald arrays
  // Set to the maximum value
  cub::DeviceReduce::Sum(vars->cub_reduce_storage,
                         vars->cub_reduce_storage_size, vars->gpu_recipEnergies,
                         vars->gpu_finalVal, imageTotal);
  CUMALLOC(&vars->cub_reduce_storage, vars->cub_reduce_storage_size);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

void CopyRefToNewCUDA(VariablesCUDA *vars, uint box, uint imageTotal) {
  cudaMemcpy(vars->gpu_sumRnew[box], vars->gpu_sumRref[box],
             imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
  cudaMemcpy(vars->gpu_sumInew[box], vars->gpu_sumIref[box],
             imageTotal * sizeof(double), cudaMemcpyDeviceToDevice);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

void UpdateEnergyVecs(VariablesCUDA *vars, int newVecLen, bool electrostatic) {
  // If we haven't exceeded the previous maximum size, we can reuse the storage
  if (vars->gpu_energyVecLen >= newVecLen)
    return;

  // Free the current allocations if this isn't the first allocation
  if (vars->gpu_energyVecLen > 0) {
    CUFREE(vars->gpu_LJEn);
    CUFREE(vars->gpu_vT11);
    CUFREE(vars->gpu_vT12);
    CUFREE(vars->gpu_vT13);
    CUFREE(vars->gpu_vT22);
    CUFREE(vars->gpu_vT23);
    CUFREE(vars->gpu_vT33);
    if (electrostatic) {
      CUFREE(vars->gpu_REn);
      CUFREE(vars->gpu_rT11);
      CUFREE(vars->gpu_rT12);
      CUFREE(vars->gpu_rT13);
      CUFREE(vars->gpu_rT22);
      CUFREE(vars->gpu_rT23);
      CUFREE(vars->gpu_rT33);
    }
  }
  vars->gpu_energyVecLen = newVecLen;
  CUMALLOC((void **)&vars->gpu_LJEn, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT11, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT12, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT13, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT22, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT23, newVecLen * sizeof(double));
  CUMALLOC((void **)&vars->gpu_vT33, newVecLen * sizeof(double));
  if (electrostatic) {
    CUMALLOC((void **)&vars->gpu_REn, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT11, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT12, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT13, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT22, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT23, newVecLen * sizeof(double));
    CUMALLOC((void **)&vars->gpu_rT33, newVecLen * sizeof(double));
  }

  // Check if more temporary storage is needed for this larger reduction size.
  // If so, free and malloc a new array for the additional space.
  void *d_temp_storage = nullptr;
  size_t temp_storage_bytes = 0;
  cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_LJEn,
                         vars->gpu_finalVal, vars->gpu_energyVecLen);
  if (temp_storage_bytes > vars->cub_energyVec_storage_size) {
    // Free the current allocation if this isn't the first allocation
    if (vars->cub_energyVec_storage_size > 0) {
      CUFREE(vars->cub_energyVec_storage);
    }
    vars->cub_energyVec_storage_size = temp_storage_bytes;
    CUMALLOC(&(vars->cub_energyVec_storage), vars->cub_energyVec_storage_size);
  }
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
  CUFREE(vars->gpu_wT11);
  CUFREE(vars->gpu_wT12);
  CUFREE(vars->gpu_wT13);
  CUFREE(vars->gpu_wT22);
  CUFREE(vars->gpu_wT23);
  CUFREE(vars->gpu_wT33);
  CUFREE(vars->gpu_recipEnergies);
  CUFREE(vars->gpu_particleUsed);
  CUFREE(vars->cub_reduce_storage);

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

void DestroyExp6CUDAVars(VariablesCUDA *vars) {
  CUFREE(vars->gpu_rMin);
  CUFREE(vars->gpu_rMaxSq);
  CUFREE(vars->gpu_expConst);
}

void DestroyCUDAVars(VariablesCUDA *vars) {
  CUFREE(vars->cub_energyVec_storage);
  CUFREE(vars->gpu_sigmaSq);
  CUFREE(vars->gpu_epsilon_Cn);
  CUFREE(vars->gpu_n);
  CUFREE(vars->gpu_VDW_Kind);
  CUFREE(vars->gpu_isMartini);
  CUFREE(vars->gpu_count);
  CUFREE(vars->gpu_rCut);
  CUFREE(vars->gpu_rCutSq);
  CUFREE(vars->gpu_rCutCoulomb);
  CUFREE(vars->gpu_rCutCoulombSq);
  CUFREE(vars->gpu_rCutLow);
  CUFREE(vars->gpu_rOn);
  CUFREE(vars->gpu_alpha);
  CUFREE(vars->gpu_alphaSq);
  CUFREE(vars->gpu_ewald);
  CUFREE(vars->gpu_diElectric_1);
  CUFREE(vars->gpu_finalVal);
  CUFREE(vars->gpu_x);
  CUFREE(vars->gpu_y);
  CUFREE(vars->gpu_z);
  CUFREE(vars->gpu_molCharge);
  CUFREE(vars->gpu_particleCharge);
  CUFREE(vars->gpu_particleKind);
  CUFREE(vars->gpu_particleMol);
  CUFREE(vars->gpu_dx);
  CUFREE(vars->gpu_dy);
  CUFREE(vars->gpu_dz);
  CUFREE(vars->gpu_nx);
  CUFREE(vars->gpu_ny);
  CUFREE(vars->gpu_nz);
  CUFREE(vars->gpu_comx);
  CUFREE(vars->gpu_comy);
  CUFREE(vars->gpu_comz);
  CUFREE(vars->gpu_LJEn);
  CUFREE(vars->gpu_REn);
  CUFREE(vars->gpu_rT11);
  CUFREE(vars->gpu_rT12);
  CUFREE(vars->gpu_rT13);
  CUFREE(vars->gpu_rT22);
  CUFREE(vars->gpu_rT23);
  CUFREE(vars->gpu_rT33);
  CUFREE(vars->gpu_vT11);
  CUFREE(vars->gpu_vT12);
  CUFREE(vars->gpu_vT13);
  CUFREE(vars->gpu_vT22);
  CUFREE(vars->gpu_vT23);
  CUFREE(vars->gpu_vT33);
  CUFREE(vars->gpu_rt_k_x);
  CUFREE(vars->gpu_rt_k_y);
  CUFREE(vars->gpu_rt_k_z);
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

  // free GPU memory for lambda variables
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
