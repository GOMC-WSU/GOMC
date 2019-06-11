/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA

#include <cuda.h>
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "cub/cub.cuh"
#include <stdio.h>

using namespace cub;

void CallBoxInterForceGPU(VariablesCUDA *vars,
                          vector<uint> &pair1,
                          vector<uint> &pair2,
                          XYZArray const &currentCoords,
                          XYZArray const &currentCOM,
                          BoxDimensions const &boxAxes,
                          bool electrostatic,
                          vector<real> &particleCharge,
                          vector<int> &particleKind,
                          vector<int> &particleMol,
                          real &rT11,
                          real &rT12,
                          real &rT13,
                          real &rT22,
                          real &rT23,
                          real &rT33,
                          real &vT11,
                          real &vT12,
                          real &vT13,
                          real &vT22,
                          real &vT23,
                          real &vT33,
                          uint const box)
{
  int atomNumber = currentCoords.Count();
  int molNumber = currentCOM.Count();
  int *gpu_pair1, *gpu_pair2;
  int *gpu_particleKind;
  int *gpu_particleMol;
  int blocksPerGrid, threadsPerBlock;
  real *gpu_particleCharge;
  real *gpu_final_value;

  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_particleKind, particleKind.size() * sizeof(int));
  cudaMalloc((void**) &gpu_particleMol, particleMol.size() * sizeof(int));
  cudaMalloc((void**) &gpu_final_value, sizeof(real));

  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, currentCOM.x, molNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, currentCOM.y, molNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, currentCOM.z, molNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0],
             particleKind.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0],
             particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);

  // Run the kernel...
  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size() / threadsPerBlock) + 1;
  BoxInterForceGPU <<< blocksPerGrid, threadsPerBlock>>>(gpu_pair1,
      gpu_pair2,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      vars->gpu_comx,
      vars->gpu_comy,
      vars->gpu_comz,
      boxAxes.GetAxis(box).x,
      boxAxes.GetAxis(box).y,
      boxAxes.GetAxis(box).z,
      electrostatic,
      gpu_particleCharge,
      gpu_particleKind,
      gpu_particleMol,
      vars->gpu_rT11,
      vars->gpu_rT12,
      vars->gpu_rT13,
      vars->gpu_rT22,
      vars->gpu_rT23,
      vars->gpu_rT33,
      vars->gpu_vT11,
      vars->gpu_vT12,
      vars->gpu_vT13,
      vars->gpu_vT22,
      vars->gpu_vT23,
      vars->gpu_vT33,
      pair1.size(),
      vars->gpu_sigmaSq,
      vars->gpu_epsilon_Cn,
      vars->gpu_n,
      vars->gpu_VDW_Kind,
      vars->gpu_isMartini,
      vars->gpu_count,
      vars->gpu_rCut,
      vars->gpu_rCutCoulomb,
      vars->gpu_rCutLow,
      vars->gpu_rOn,
      vars->gpu_alpha,
      vars->gpu_ewald,
      vars->gpu_diElectric_1,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      vars->gpu_nonOrth,
      box);
  cudaDeviceSynchronize();
  // ReduceSum // Virial of LJ
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT11,
                    gpu_final_value, pair1.size());
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT11,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT11, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT12,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT12, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT13,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT13, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT22,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT22, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT23,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT23, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_vT33,
                    gpu_final_value, pair1.size());
  cudaMemcpy(&vT33, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);

  if(electrostatic) {
    // ReduceSum // Virial of Coulomb
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT11,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT11, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT12,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT12, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT13,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT13, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT22,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT22, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT23,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT23, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT33,
                      gpu_final_value, pair1.size());
    cudaMemcpy(&rT33, gpu_final_value, sizeof(real),
               cudaMemcpyDeviceToHost);
  }

  cudaFree(d_temp_storage);
  cudaFree(gpu_pair1);
  cudaFree(gpu_pair2);
  cudaFree(gpu_particleKind);
  cudaFree(gpu_particleMol);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_final_value);
}

void CallForceReciprocalGPU(VariablesCUDA *vars,
                            XYZArray const &currentCoords,
                            XYZArray const &currentCOMDiff,
                            vector<real> &particleCharge,
                            real &rT11,
                            real &rT12,
                            real &rT13,
                            real &rT22,
                            real &rT23,
                            real &rT33,
                            uint imageSize,
                            real constVal,
                            uint box)
{
  int atomNumber = currentCoords.Count();
  int blocksPerGrid, threadsPerBlock;
  real *gpu_particleCharge;
  real *gpu_final_value;

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_final_value, sizeof(real));

  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dx, currentCOMDiff.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dy, currentCOMDiff.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dz, currentCOMDiff.z, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);

  // Run the kernel...
  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  ForceReciprocalGPU <<< blocksPerGrid,
                     threadsPerBlock>>>(vars->gpu_x,
                                        vars->gpu_y,
                                        vars->gpu_z,
                                        vars->gpu_dx,
                                        vars->gpu_dy,
                                        vars->gpu_dz,
                                        vars->gpu_kxRef[box],
                                        vars->gpu_kyRef[box],
                                        vars->gpu_kzRef[box],
                                        vars->gpu_prefactRef[box],
                                        vars->gpu_hsqrRef[box],
                                        vars->gpu_sumRref[box],
                                        vars->gpu_sumIref[box],
                                        gpu_particleCharge,
                                        vars->gpu_rT11,
                                        vars->gpu_rT12,
                                        vars->gpu_rT13,
                                        vars->gpu_rT22,
                                        vars->gpu_rT23,
                                        vars->gpu_rT33,
                                        constVal,
                                        imageSize,
                                        atomNumber);

  // ReduceSum // Virial of Reciprocal
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT11,
                    gpu_final_value, imageSize);
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT11,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT11, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT12,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT12, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT13,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT13, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT22,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT22, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT23,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT23, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, vars->gpu_rT33,
                    gpu_final_value, imageSize);
  cudaMemcpy(&rT33, gpu_final_value, sizeof(real),
             cudaMemcpyDeviceToHost);

  cudaFree(gpu_particleCharge);
  cudaFree(gpu_final_value);
  cudaFree(d_temp_storage);
}

__global__ void BoxInterForceGPU(int *gpu_pair1,
                                 int *gpu_pair2,
                                 real *gpu_x,
                                 real *gpu_y,
                                 real *gpu_z,
                                 real *gpu_comx,
                                 real *gpu_comy,
                                 real *gpu_comz,
                                 real xAxes,
                                 real yAxes,
                                 real zAxes,
                                 bool electrostatic,
                                 real *gpu_particleCharge,
                                 int *gpu_particleKind,
                                 int *gpu_particleMol,
                                 real *gpu_rT11,
                                 real *gpu_rT12,
                                 real *gpu_rT13,
                                 real *gpu_rT22,
                                 real *gpu_rT23,
                                 real *gpu_rT33,
                                 real *gpu_vT11,
                                 real *gpu_vT12,
                                 real *gpu_vT13,
                                 real *gpu_vT22,
                                 real *gpu_vT23,
                                 real *gpu_vT33,
                                 int pairSize,
                                 real *gpu_sigmaSq,
                                 real *gpu_epsilon_Cn,
                                 real *gpu_n,
                                 int *gpu_VDW_Kind,
                                 int *gpu_isMartini,
                                 int *gpu_count,
                                 real *gpu_rCut,
                                 real *gpu_rCutCoulomb,
                                 real *gpu_rCutLow,
                                 real *gpu_rOn,
                                 real *gpu_alpha,
                                 int *gpu_ewald,
                                 real *gpu_diElectric_1,
                                 real *gpu_cell_x,
                                 real *gpu_cell_y,
                                 real *gpu_cell_z,
                                 real *gpu_Invcell_x,
                                 real *gpu_Invcell_y,
                                 real *gpu_Invcell_z,
                                 int *gpu_nonOrth,
                                 int box)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= pairSize)
    return;

  real distSq;
  real virX, virY, virZ;
  real pRF = 0.0, qi_qj, pVF = 0.0;
  //tensors for VDW and real part of electrostatic
  gpu_vT11[threadID] = 0.0, gpu_vT22[threadID] = 0.0, gpu_vT33[threadID] = 0.0;
  gpu_rT11[threadID] = 0.0, gpu_rT22[threadID] = 0.0, gpu_rT33[threadID] = 0.0;
  // extra tensors reserved for later on
  gpu_vT12[threadID] = 0.0, gpu_vT13[threadID] = 0.0, gpu_vT23[threadID] = 0.0;
  gpu_rT12[threadID] = 0.0, gpu_rT13[threadID] = 0.0, gpu_rT23[threadID] = 0.0;
  real diff_comx, diff_comy, diff_comz;
  real cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);

  if(InRcutGPU(distSq, virX, virY, virZ, gpu_x[gpu_pair1[threadID]],
               gpu_y[gpu_pair1[threadID]], gpu_z[gpu_pair1[threadID]],
               gpu_x[gpu_pair2[threadID]], gpu_y[gpu_pair2[threadID]],
               gpu_z[gpu_pair2[threadID]], xAxes, yAxes, zAxes, xAxes / 2.0,
               yAxes / 2.0, zAxes / 2.0, cutoff, gpu_nonOrth[0],
               gpu_cell_x, gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
               gpu_Invcell_z)) {
    diff_comx = gpu_comx[gpu_particleMol[gpu_pair1[threadID]]] -
                gpu_comx[gpu_particleMol[gpu_pair2[threadID]]];
    diff_comy = gpu_comy[gpu_particleMol[gpu_pair1[threadID]]] -
                gpu_comy[gpu_particleMol[gpu_pair2[threadID]]];
    diff_comz = gpu_comz[gpu_particleMol[gpu_pair1[threadID]]] -
                gpu_comz[gpu_particleMol[gpu_pair2[threadID]]];

    diff_comx = MinImageSignedGPU(diff_comx, xAxes, xAxes / 2.0);
    diff_comy = MinImageSignedGPU(diff_comy, yAxes, yAxes / 2.0);
    diff_comz = MinImageSignedGPU(diff_comz, zAxes, zAxes / 2.0);

    if(electrostatic) {
      qi_qj = gpu_particleCharge[gpu_pair1[threadID]] *
              gpu_particleCharge[gpu_pair2[threadID]];
      pRF = CalcCoulombForceGPU(distSq, qi_qj, gpu_VDW_Kind[0], gpu_ewald[0],
                                gpu_isMartini[0], gpu_alpha[box], gpu_rCutCoulomb[box],
                                gpu_diElectric_1[0]);

      gpu_rT11[threadID] = pRF * (virX * diff_comx);
      gpu_rT22[threadID] = pRF * (virY * diff_comy);
      gpu_rT33[threadID] = pRF * (virZ * diff_comz);

      //extra tensor calculations
      gpu_rT12[threadID] = pRF * (0.5 * (virX * diff_comy + virY * diff_comx));
      gpu_rT13[threadID] = pRF * (0.5 * (virX * diff_comz + virZ * diff_comx));
      gpu_rT23[threadID] = pRF * (0.5 * (virY * diff_comz + virZ * diff_comy));
    }

    pVF = CalcEnForceGPU(distSq, gpu_particleKind[gpu_pair1[threadID]],
                         gpu_particleKind[gpu_pair2[threadID]],
                         gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_rCut[0],
                         gpu_rOn[0], gpu_isMartini[0], gpu_VDW_Kind[0],
                         gpu_count[0]);

    gpu_vT11[threadID] = pVF * (virX * diff_comx);
    gpu_vT22[threadID] = pVF * (virY * diff_comy);
    gpu_vT33[threadID] = pVF * (virZ * diff_comz);

    //extra tensor calculations
    gpu_vT12[threadID] = pVF * (0.5 * (virX * diff_comy + virY * diff_comx));
    gpu_vT13[threadID] = pVF * (0.5 * (virX * diff_comz + virZ * diff_comx));
    gpu_vT23[threadID] = pVF * (0.5 * (virY * diff_comz + virZ * diff_comy));
  }
}


__global__ void ForceReciprocalGPU(real *gpu_x,
                                   real *gpu_y,
                                   real *gpu_z,
                                   real *gpu_comDx,
                                   real *gpu_comDy,
                                   real *gpu_comDz,
                                   real *gpu_kxRef,
                                   real *gpu_kyRef,
                                   real *gpu_kzRef,
                                   real *gpu_prefactRef,
                                   real *gpu_hsqrRef,
                                   real *gpu_sumRref,
                                   real *gpu_sumIref,
                                   real *gpu_particleCharge,
                                   real *gpu_rT11,
                                   real *gpu_rT12,
                                   real *gpu_rT13,
                                   real *gpu_rT22,
                                   real *gpu_rT23,
                                   real *gpu_rT33,
                                   real constVal,
                                   uint imageSize,
                                   uint atomNumber)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  real factor, arg;
  int i;
  factor = gpu_prefactRef[threadID] * (gpu_sumRref[threadID] *
                                       gpu_sumRref[threadID] +
                                       gpu_sumIref[threadID] *
                                       gpu_sumIref[threadID]);
  gpu_rT11[threadID] = factor * (1.0 - 2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kxRef[threadID] * gpu_kxRef[threadID]);
  gpu_rT12[threadID] = factor * (-2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kxRef[threadID] * gpu_kyRef[threadID]);
  gpu_rT13[threadID] = factor * (-2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kxRef[threadID] * gpu_kzRef[threadID]);
  gpu_rT22[threadID] = factor * (1.0 - 2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kyRef[threadID] * gpu_kyRef[threadID]);
  gpu_rT23[threadID] = factor * (-2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kyRef[threadID] * gpu_kzRef[threadID]);
  gpu_rT33[threadID] = factor * (1.0 - 2.0 *
                                 (constVal + 1.0 / gpu_hsqrRef[threadID]) *
                                 gpu_kzRef[threadID] * gpu_kzRef[threadID]);

  //Intramolecular part
  for(i = 0; i < atomNumber; i++) {
    arg = DotProductGPU(gpu_kxRef[threadID], gpu_kyRef[threadID],
                        gpu_kzRef[threadID], gpu_x[i], gpu_y[i], gpu_z[i]);

    factor = gpu_prefactRef[threadID] * 2.0 *
             (gpu_sumIref[threadID] * cos(arg) - gpu_sumRref[threadID] * sin(arg)) *
             gpu_particleCharge[i];

    gpu_rT11[threadID] += factor * (gpu_kxRef[threadID] * gpu_comDx[i]);
    gpu_rT12[threadID] += factor * 0.5 * (gpu_kxRef[threadID] * gpu_comDy[i] +
                                          gpu_kyRef[threadID] * gpu_comDx[i]);
    gpu_rT13[threadID] += factor * 0.5 * (gpu_kxRef[threadID] * gpu_comDz[i] +
                                          gpu_kzRef[threadID] * gpu_comDx[i]);
    gpu_rT22[threadID] += factor * (gpu_kyRef[threadID] * gpu_comDy[i]);
    gpu_rT13[threadID] += factor * 0.5 * (gpu_kyRef[threadID] * gpu_comDz[i] +
                                          gpu_kzRef[threadID] * gpu_comDy[i]);
    gpu_rT33[threadID] += factor * (gpu_kzRef[threadID] * gpu_comDz[i]);
  }
}

__device__ real CalcCoulombForceGPU(real distSq, real qi_qj,
                                      int gpu_VDW_Kind, int gpu_ewald,
                                      int gpu_isMartini, real gpu_alpha,
                                      real gpu_rCutCoulomb, real gpu_diElectric_1)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombVirParticleGPU(distSq, qi_qj, gpu_alpha);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcCoulombVirShiftGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                          gpu_rCutCoulomb, gpu_diElectric_1);
  } else
    return CalcCoulombVirSwitchGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                   gpu_rCutCoulomb);
}

__device__ real CalcEnForceGPU(real distSq, int kind1, int kind2,
                                 real *gpu_sigmaSq, real *gpu_n,
                                 real *gpu_epsilon_Cn, real gpu_rCut,
                                 real gpu_rOn, int gpu_isMartini,
                                 int gpu_VDW_Kind, int gpu_count)
{
  if((gpu_rCut * gpu_rCut) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcVirParticleGPU(distSq, index, gpu_sigmaSq, gpu_n,
                              gpu_epsilon_Cn);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcVirShiftGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcVirSwitchMartiniGPU(distSq, index, gpu_sigmaSq, gpu_n,
                                   gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  } else
    return CalcVirSwitchGPU(distSq, index, gpu_sigmaSq, gpu_epsilon_Cn, gpu_n,
                            gpu_rCut, gpu_rOn);
}

//ElectroStatic Calculation
//**************************************************************//
__device__ real CalcCoulombVirParticleGPU(real distSq, real qi_qj,
    real gpu_alpha)
{
  real dist = sqrt(distSq);
  real constValue = 2.0 * gpu_alpha / sqrt(M_PI);
  real expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
  real temp = 1.0 - erf(gpu_alpha * dist);
  return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
}

__device__ real CalcCoulombVirShiftGPU(real distSq, real qi_qj,
    int gpu_ewald, real gpu_alpha)
{
  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real constValue = 2.0 * gpu_alpha / sqrt(M_PI);
    real expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    real temp = 1.0 - erf(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    real dist = sqrt(distSq);
    return qi_qj / (distSq * dist);
  }
}
__device__ real CalcCoulombVirSwitchMartiniGPU(real distSq, real qi_qj,
    int gpu_ewald,
    real gpu_alpha,
    real gpu_rCut,
    real gpu_diElectric_1)
{
  real power2 = 2.0;
  real power3 = 3.0;

  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real constValue = 2.0 * gpu_alpha / sqrt(M_PI);
    real expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    real temp = 1.0 - erf(gpu_alpha * dist);
    return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    real dist = sqrt(distSq);
    real rij_ronCoul_2 = distSq;
    real rij_ronCoul_3 = dist * distSq;

    real A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, power3) *
                pow(gpu_rCut, power2));
    real B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, power3) *
                pow(gpu_rCut, power3));

    real virCoul = A1 / rij_ronCoul_2 + B1 / rij_ronCoul_3;
    return qi_qj * gpu_diElectric_1 * ( 1.0 / (dist * distSq) + virCoul / dist);
  }
}

__device__ real CalcCoulombVirSwitchGPU(real distSq, real qi_qj,
    int gpu_ewald, real gpu_alpha,
    real gpu_rCut)
{
  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real constValue = 2.0 * gpu_alpha / sqrt(M_PI);
    real expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    real temp = 1.0 - erf(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    real rCutSq = gpu_rCut * gpu_rCut;
    real dist = sqrt(distSq);
    real switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;

    real dSwitchVal = 2.0 * (distSq / rCutSq - 1.0) * 2.0 * dist / rCutSq;
    return -1.0 * qi_qj * (dSwitchVal / distSq - switchVal / (distSq * dist));
  }
}

//VDW Calculation
//*****************************************************************//
__device__ real CalcVirParticleGPU(real distSq, int index,
                                     real *gpu_sigmaSq, real *gpu_n,
                                     real *gpu_epsilon_Cn)
{
  real power2 = 2.0;

  real rNeg2 = 1.0 / distSq;
  real rRat2 = gpu_sigmaSq[index] * rNeg2;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
  real repulse = pow(rRat2, gpu_n[index] / power2);
  return gpu_epsilon_Cn[index] * 6.0 *
         ((gpu_n[index] / 6.0) * repulse - attract) * rNeg2;
}

__device__ real CalcVirShiftGPU(real distSq, int index, real *gpu_sigmaSq,
                                  real *gpu_n, real *gpu_epsilon_Cn)
{
  real rNeg2 = 1.0 / distSq;
  real rRat2 = gpu_sigmaSq[index] * rNeg2;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
  real repulse = pow(rRat2, gpu_n[index] / 2.0);
  return gpu_epsilon_Cn[index] * 6.0 *
         ((gpu_n[index] / 6.0) * repulse - attract) * rNeg2;
}

__device__ real CalcVirSwitchMartiniGPU(real distSq, int index,
    real *gpu_sigmaSq, real *gpu_n,
    real *gpu_epsilon_Cn,
    real gpu_rCut, real gpu_rOn)
{
  real power2 = 2.0;
  real power3 = 3.0;
  real power8 = 8.0;

  real r_1 = 1.0 / sqrt(distSq);
  real r_8 = pow(r_1, power8);
  real r_n2 = pow(r_1, gpu_n[index] + 2);

  real rij_ron = sqrt(distSq) - gpu_rOn;
  real rij_ron_2 = rij_ron * rij_ron;
  real rij_ron_3 = rij_ron_2 * rij_ron;

  real pn = gpu_n[index];
  real An = pn * ((pn + 1) * gpu_rOn - (pn + 4) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, power2));
  real Bn = -pn * ((pn + 1) * gpu_rOn - (pn + 3) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, power3));

  real sig6 = pow(gpu_sigmaSq[index], power3);
  real sign = pow(gpu_sigmaSq[index], pn / 2);

  real A6 = 6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 4) * gpu_rCut) /
              (pow(gpu_rCut, power8) * pow(gpu_rCut - gpu_rOn, power2));
  real B6 = -6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 3) * gpu_rCut) /
              (pow(gpu_rCut, power8) * pow(gpu_rCut - gpu_rOn, power3));

  real dshifttempRep = An * rij_ron_2 + Bn * rij_ron_3;
  real dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;

  const real dshiftRep = ( distSq > gpu_rOn * gpu_rOn ?
                             dshifttempRep * r_1 : 0);
  const real dshiftAtt = ( distSq > gpu_rOn * gpu_rOn ?
                             dshifttempAtt * r_1 : 0);
  real Wij = gpu_epsilon_Cn[index] * (sign * (pn * r_n2 + dshiftRep) -
                                        sig6 * (6.0 * r_8 + dshiftAtt));
  return Wij;
}

__device__ real CalcVirSwitchGPU(real distSq, int index,
                                   real *gpu_sigmaSq, real *gpu_epsilon_Cn,
                                   real *gpu_n, real gpu_rCut,
                                   real gpu_rOn)
{
  real powerneg3 = -3.0;

  real rCutSq = gpu_rCut * gpu_rCut;
  real rCutSq_rijSq = rCutSq - distSq;
  real rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;
  real rOnSq = gpu_rOn * gpu_rOn;

  real rNeg2 = 1.0 / distSq;
  real rRat2 = rNeg2 * gpu_sigmaSq[index];
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
  real repulse = pow(rRat2, gpu_n[index] / 2.0);
  real factor1 = rCutSq - 3 * rOnSq;
  real factor2 = pow((rCutSq - rOnSq), powerneg3);

  real fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
  real fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

  const real factE = ( distSq > rOnSq ? fE : 1.0);
  const real factW = ( distSq > rOnSq ? fW : 0.0);

  real Wij = gpu_epsilon_Cn[index] * 6.0 *
               ((gpu_n[index] / 6.0) * repulse - attract) * rNeg2;
  real Eij = gpu_epsilon_Cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

#endif
