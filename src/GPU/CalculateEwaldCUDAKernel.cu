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

#include <vector>

#include "BoxDimensions.h"
#include "CUDAMemoryManager.cuh"
#include "CalculateEwaldCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "cub/cub.cuh"

using namespace cub;

#define IMAGES_PER_BLOCK 64
#define PARTICLES_PER_BLOCK 64

#define FULL_MASK 0xffffffff

// Use this function when calculating the reciprocal terms
// for a new volume. A change in box dimensions.
void CallBoxReciprocalSetupGPU(VariablesCUDA *vars, XYZArray const &coords,
                               double const *kx, double const *ky,
                               double const *kz,
                               const std::vector<double> &particleCharge,
                               uint imageSize, double *sumRnew, double *sumInew,
                               double *prefact, double *hsqr,
                               double &energyRecip, uint box) {
  double *gpu_particleCharge;
  double *gpu_energyRecip;
  double *gpu_final_energyRecip;
  int atomNumber = coords.Count();

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_energyRecip, imageSize * sizeof(double));
  CUMALLOC((void **)&gpu_final_energyRecip, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_kx[box], kx, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_ky[box], ky, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_kz[box], kz, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_prefact[box], prefact, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_hsqr[box], hsqr, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemset(vars->gpu_sumRnew[box], 0, imageSize * sizeof(double));
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemset(vars->gpu_sumInew[box], 0, imageSize * sizeof(double));
  checkLastErrorCUDA(__FILE__, __LINE__);

  dim3 threadsPerBlock(256, 1, 1);
  dim3 blocksPerGrid((int)(imageSize / threadsPerBlock.x) + 1,
                     (int)(atomNumber / PARTICLES_PER_BLOCK) + 1, 1);
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kx[box],
      vars->gpu_ky[box], vars->gpu_kz[box], atomNumber, gpu_particleCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  // Need just one thread per image for this kernel.
  blocksPerGrid.y = 1;
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_prefact[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      gpu_energyRecip, imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  cudaMemcpy(&energyRecip, gpu_final_energyRecip, sizeof(double),
             cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecip);
  CUFREE(gpu_final_energyRecip);
  CUFREE(d_temp_storage);
}

// Use this function when calculating the reciprocal terms
// with the current volume.
void CallBoxReciprocalSumsGPU(VariablesCUDA *vars, XYZArray const &coords,
                              const std::vector<double> &particleCharge,
                              uint imageSize, double *sumRnew, double *sumInew,
                              double &energyRecip, uint box) {
  double *gpu_particleCharge;
  double *gpu_energyRecip;
  double *gpu_final_energyRecip;
  int atomNumber = coords.Count();

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_energyRecip, imageSize * sizeof(double));
  CUMALLOC((void **)&gpu_final_energyRecip, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemset(vars->gpu_sumRnew[box], 0, imageSize * sizeof(double));
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemset(vars->gpu_sumInew[box], 0, imageSize * sizeof(double));
  checkLastErrorCUDA(__FILE__, __LINE__);

  dim3 threadsPerBlock(256, 1, 1);
  dim3 blocksPerGrid((int)(imageSize / threadsPerBlock.x) + 1,
                     (int)(atomNumber / PARTICLES_PER_BLOCK) + 1, 1);
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], atomNumber,
      gpu_particleCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  // Need just one thread per image for this kernel.
  blocksPerGrid.y = 1;
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_prefactRef[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      gpu_energyRecip, imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  cudaMemcpy(&energyRecip, gpu_final_energyRecip, sizeof(double),
             cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecip);
  CUFREE(gpu_final_energyRecip);
  CUFREE(d_temp_storage);
}

__global__ void BoxReciprocalSumsGPU(double *gpu_x, double *gpu_y,
                                     double *gpu_z, double *gpu_kx,
                                     double *gpu_ky, double *gpu_kz,
                                     int atomNumber, double *gpu_particleCharge,
                                     double *gpu_sumRnew, double *gpu_sumInew,
                                     int imageSize) {
  __shared__ double shared_coords[PARTICLES_PER_BLOCK * 3];
  int imageID = blockIdx.x * blockDim.x + threadIdx.x;
  int offset_coordinates_index = blockIdx.y * PARTICLES_PER_BLOCK;
  int numberOfAtoms =
      min(PARTICLES_PER_BLOCK, atomNumber - offset_coordinates_index);
  double sumR = 0.0, sumI = 0.0;

  if (threadIdx.x < numberOfAtoms) {
    shared_coords[threadIdx.x * 3] =
        gpu_x[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 3 + 1] =
        gpu_y[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 3 + 2] =
        gpu_z[offset_coordinates_index + threadIdx.x];
  }

  if (imageID >= imageSize)
    return;

  __syncthreads();
  for (int particleID = 0; particleID < numberOfAtoms; particleID++) {
    double dot = DotProductGPU(gpu_kx[imageID], gpu_ky[imageID],
                               gpu_kz[imageID], shared_coords[particleID * 3],
                               shared_coords[particleID * 3 + 1],
                               shared_coords[particleID * 3 + 2]);
    double dotsin, dotcos;
    sincos(dot, &dotsin, &dotcos);
    sumR += gpu_particleCharge[offset_coordinates_index + particleID] * dotcos;
    sumI += gpu_particleCharge[offset_coordinates_index + particleID] * dotsin;
  }

  atomicAdd(&gpu_sumRnew[imageID], sumR);
  atomicAdd(&gpu_sumInew[imageID], sumI);
}

__global__ void BoxReciprocalGPU(double *gpu_prefact, double *gpu_sumRnew,
                                 double *gpu_sumInew, double *gpu_energyRecip,
                                 int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  gpu_energyRecip[threadID] = ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
                                gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
                               gpu_prefact[threadID]);
}

void CallMolReciprocalGPU(VariablesCUDA *vars, XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          const std::vector<double> &particleCharge,
                          uint imageSize, double *sumRnew, double *sumInew,
                          double &energyRecipNew, uint box) {
  // Calculate atom number
  int atomNumber = currentCoords.Count();
  int newCoordsNumber = newCoords.Count();
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_energyRecipNew, imageSize * sizeof(double));
  CUMALLOC((void **)&gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_nx, newCoords.x, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ny, newCoords.y, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nz, newCoords.z, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  MolReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_nx, vars->gpu_ny,
      vars->gpu_nz, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], atomNumber, gpu_particleCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], vars->gpu_sumRref[box],
      vars->gpu_sumIref[box], vars->gpu_prefactRef[box], gpu_energyRecipNew,
      imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew, sizeof(double),
             cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecipNew);
  CUFREE(gpu_final_energyRecipNew);
  CUFREE(d_temp_storage);
}

// Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &particleCharge,
                                      uint imageSize, double *sumRnew,
                                      double *sumInew, double &energyRecipNew,
                                      const double lambdaCoef, uint box) {
  // Calculate atom number
  int atomNumber = coords.Count();
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_energyRecipNew, imageSize * sizeof(double));
  CUMALLOC((void **)&gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  ChangeLambdaMolReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], atomNumber,
      gpu_particleCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_sumRref[box], vars->gpu_sumIref[box], vars->gpu_prefactRef[box],
      gpu_energyRecipNew, lambdaCoef, imageSize);

  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew, sizeof(double),
             cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecipNew);
  CUFREE(gpu_final_energyRecipNew);
  CUFREE(d_temp_storage);
}

void CallSwapReciprocalGPU(VariablesCUDA *vars, XYZArray const &coords,
                           const std::vector<double> &particleCharge,
                           uint imageSize, double *sumRnew, double *sumInew,
                           const bool insert, double &energyRecipNew,
                           uint box) {
  // Calculate atom number
  int atomNumber = coords.Count();
  // given coordinates
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_energyRecipNew, imageSize * sizeof(double));
  CUMALLOC((void **)&gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  SwapReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], atomNumber,
      gpu_particleCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_sumRref[box], vars->gpu_sumIref[box], vars->gpu_prefactRef[box],
      insert, gpu_energyRecipNew, imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
  //#ifndef NDEBUG
  // In the future maybe we could remove this for Nondebug?
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  //#endif

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew, sizeof(double),
             cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecipNew);
  CUFREE(gpu_final_energyRecipNew);
  CUFREE(d_temp_storage);
}

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars, uint imageSize,
                                  double *sumRnew, double *sumInew, uint box) {
  cudaMemcpy(vars->gpu_sumRnew[box], sumRnew, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_sumInew[box], sumInew, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
}

void CallBoxForceReciprocalGPU(
    VariablesCUDA *vars, XYZArray &atomForceRec, XYZArray &molForceRec,
    const std::vector<double> &particleCharge,
    const std::vector<int> &particleMol,
    const std::vector<bool> &particleHasNoCharge, const bool *particleUsed,
    const std::vector<int> &startMol, const std::vector<int> &lengthMol,
    double alpha, double alphaSq, double constValue, uint imageSize,
    XYZArray const &molCoords, BoxDimensions const &boxAxes, int box) {
  int atomCount = atomForceRec.Count();
  int molCount = molForceRec.Count();
  double *gpu_particleCharge;
  int *gpu_particleMol;
  bool *gpu_particleHasNoCharge, *gpu_particleUsed;
  bool *arr_particleHasNoCharge = new bool[particleHasNoCharge.size()];
  int *gpu_startMol, *gpu_lengthMol;

  // particleHasNoCharge is stored in vector<bool>, so in order to copy it to
  // GPU it needs to be stored in bool[]. because: std::vector<bool> : Does not
  // necessarily store its elements as a contiguous array
  for (int i = 0; i < particleHasNoCharge.size(); i++) {
    arr_particleHasNoCharge[i] = particleHasNoCharge[i];
  }

  // calculate block and grid sizes
  dim3 threadsPerBlock(256, 1, 1);
  int blocksPerGridX = (int)(atomCount / threadsPerBlock.x) + 1;
  int blocksPerGridY = (int)(imageSize / IMAGES_PER_BLOCK) + 1;
  dim3 blocksPerGrid(blocksPerGridX, blocksPerGridY, 1);

  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_particleHasNoCharge,
           particleHasNoCharge.size() * sizeof(bool));
  CUMALLOC((void **)&gpu_particleUsed, atomCount * sizeof(bool));
  CUMALLOC((void **)&gpu_startMol, startMol.size() * sizeof(int));
  CUMALLOC((void **)&gpu_lengthMol, lengthMol.size() * sizeof(int));
  CUMALLOC((void **)&gpu_particleMol, particleMol.size() * sizeof(int));

  cudaMemcpy(vars->gpu_aForceRecx, atomForceRec.x, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForceRecy, atomForceRec.y, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForceRecz, atomForceRec.z, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecx, molForceRec.x, sizeof(double) * molCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecy, molForceRec.y, sizeof(double) * molCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecz, molForceRec.z, sizeof(double) * molCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             sizeof(double) * particleCharge.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], sizeof(int) * particleMol.size(),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleHasNoCharge, arr_particleHasNoCharge,
             sizeof(bool) * particleHasNoCharge.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleUsed, particleUsed, sizeof(bool) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, molCoords.x, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, molCoords.y, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, molCoords.z, sizeof(double) * atomCount,
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_startMol, &startMol[0], sizeof(int) * startMol.size(),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_lengthMol, &lengthMol[0], sizeof(int) * lengthMol.size(),
             cudaMemcpyHostToDevice);

  checkLastErrorCUDA(__FILE__, __LINE__);
  BoxForceReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_aForceRecx, vars->gpu_aForceRecy, vars->gpu_aForceRecz,
      vars->gpu_mForceRecx, vars->gpu_mForceRecy, vars->gpu_mForceRecz,
      gpu_particleCharge, gpu_particleMol, gpu_particleHasNoCharge,
      gpu_particleUsed, gpu_startMol, gpu_lengthMol, alpha, alphaSq, constValue,
      imageSize, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], vars->gpu_x, vars->gpu_y, vars->gpu_z,
      vars->gpu_prefactRef[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_isFraction, vars->gpu_molIndex, vars->gpu_lambdaCoulomb,
      vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box], vars->gpu_nonOrth, boxAxes.GetAxis(box).x,
      boxAxes.GetAxis(box).y, boxAxes.GetAxis(box).z, box, atomCount);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(atomForceRec.x, vars->gpu_aForceRecx, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(atomForceRec.y, vars->gpu_aForceRecy, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(atomForceRec.z, vars->gpu_aForceRecz, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.x, vars->gpu_mForceRecx, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.y, vars->gpu_mForceRecy, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.z, vars->gpu_mForceRecz, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  delete[] arr_particleHasNoCharge;
  CUFREE(gpu_particleCharge);
  CUFREE(gpu_particleHasNoCharge);
  CUFREE(gpu_particleUsed);
  CUFREE(gpu_startMol);
  CUFREE(gpu_lengthMol);
  CUFREE(gpu_particleMol);
}

__global__ void BoxForceReciprocalGPU(
    double *gpu_aForceRecx, double *gpu_aForceRecy, double *gpu_aForceRecz,
    double *gpu_mForceRecx, double *gpu_mForceRecy, double *gpu_mForceRecz,
    double *gpu_particleCharge, int *gpu_particleMol,
    bool *gpu_particleHasNoCharge, bool *gpu_particleUsed, int *gpu_startMol,
    int *gpu_lengthMol, double alpha, double alphaSq, double constValue,
    int imageSize, double *gpu_kx, double *gpu_ky, double *gpu_kz,
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_prefact,
    double *gpu_sumRnew, double *gpu_sumInew, bool *gpu_isFraction,
    int *gpu_molIndex, double *gpu_lambdaCoulomb, double *gpu_cell_x,
    double *gpu_cell_y, double *gpu_cell_z, double *gpu_Invcell_x,
    double *gpu_Invcell_y, double *gpu_Invcell_z, int *gpu_nonOrth, double axx,
    double axy, double axz, int box, int atomCount) {
  __shared__ double shared_kvector[IMAGES_PER_BLOCK * 3];
  int particleID = blockDim.x * blockIdx.x + threadIdx.x;
  int offset_vector_index = blockIdx.y * IMAGES_PER_BLOCK;
  int numberOfVectors = min(IMAGES_PER_BLOCK, imageSize - offset_vector_index);

  if (threadIdx.x < numberOfVectors) {
    shared_kvector[threadIdx.x * 3] = gpu_kx[offset_vector_index + threadIdx.x];
    shared_kvector[threadIdx.x * 3 + 1] =
        gpu_ky[offset_vector_index + threadIdx.x];
    shared_kvector[threadIdx.x * 3 + 2] =
        gpu_kz[offset_vector_index + threadIdx.x];
  }

  if (particleID >= atomCount || !gpu_particleUsed[particleID])
    return;
  double forceX = 0.0, forceY = 0.0, forceZ = 0.0;
  int moleculeID = gpu_particleMol[particleID];

  if (gpu_particleHasNoCharge[particleID])
    return;

  double x = gpu_x[particleID];
  double y = gpu_y[particleID];
  double z = gpu_z[particleID];
  double lambdaCoef = DeviceGetLambdaCoulomb(moleculeID, box, gpu_isFraction,
                                             gpu_molIndex, gpu_lambdaCoulomb);

  __syncthreads();
  // loop over images
  for (int vectorIndex = 0; vectorIndex < numberOfVectors; vectorIndex++) {
    double dot = x * shared_kvector[vectorIndex * 3] +
                 y * shared_kvector[vectorIndex * 3 + 1] +
                 z * shared_kvector[vectorIndex * 3 + 2];
    double dotsin, dotcos;
    sincos(dot, &dotsin, &dotcos);
    double factor = 2.0 * gpu_particleCharge[particleID] *
                    gpu_prefact[offset_vector_index + vectorIndex] *
                    lambdaCoef *
                    (dotsin * gpu_sumRnew[offset_vector_index + vectorIndex] -
                     dotcos * gpu_sumInew[offset_vector_index + vectorIndex]);

    forceX += factor * shared_kvector[vectorIndex * 3];
    forceY += factor * shared_kvector[vectorIndex * 3 + 1];
    forceZ += factor * shared_kvector[vectorIndex * 3 + 2];
  }

  // loop over other particles within the same molecule
  if (blockIdx.y == 0) {
    double intraForce = 0.0, distSq = 0.0, dist = 0.0;
    double3 distVect;
    int lastParticleWithinSameMolecule =
        gpu_startMol[particleID] + gpu_lengthMol[particleID];
    for (int otherParticle = gpu_startMol[particleID];
         otherParticle < lastParticleWithinSameMolecule; otherParticle++) {
      if (particleID != otherParticle) {
        DeviceInRcut(distSq, distVect, gpu_x, gpu_y, gpu_z, particleID,
                     otherParticle, axx, axy, axz, *gpu_nonOrth, gpu_cell_x,
                     gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                     gpu_Invcell_z);
        dist = sqrt(distSq);

        double expConstValue = exp(-1.0 * alphaSq * distSq);
        double qiqj = gpu_particleCharge[particleID] *
                      gpu_particleCharge[otherParticle] * qqFactGPU;
        intraForce = qiqj * lambdaCoef * lambdaCoef / distSq;
        intraForce *= ((erf(alpha * dist) / dist) - constValue * expConstValue);
        forceX -= intraForce * distVect.x;
        forceY -= intraForce * distVect.y;
        forceZ -= intraForce * distVect.z;
      }
    }
  }

  atomicAdd(&gpu_aForceRecx[particleID], forceX);
  atomicAdd(&gpu_aForceRecy[particleID], forceY);
  atomicAdd(&gpu_aForceRecz[particleID], forceZ);
  atomicAdd(&gpu_mForceRecx[moleculeID], forceX);
  atomicAdd(&gpu_mForceRecy[moleculeID], forceY);
  atomicAdd(&gpu_mForceRecz[moleculeID], forceZ);
}

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                  double *gpu_kx, double *gpu_ky,
                                  double *gpu_kz, int atomNumber,
                                  double *gpu_particleCharge,
                                  double *gpu_sumRnew, double *gpu_sumInew,
                                  double *gpu_sumRref, double *gpu_sumIref,
                                  double *gpu_prefactRef, const bool insert,
                                  double *gpu_energyRecipNew, int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  double sumReal = 0.0, sumImaginary = 0.0;

  for (int p = 0; p < atomNumber; p++) {
    double dotProduct =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_x[p], gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    sumReal += (gpu_particleCharge[p] * dotcos);
    sumImaginary += (gpu_particleCharge[p] * dotsin);
  }

  // If we insert the molecule to the box, we add the sum value.
  // Otherwise, we subtract the sum value
  if (insert) {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] + sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] + sumImaginary;
  } else {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] - sumImaginary;
  }

  gpu_energyRecipNew[threadID] =
      ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
        gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
       gpu_prefactRef[threadID]);
}

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber, double *gpu_particleCharge,
                                 double *gpu_sumRnew, double *gpu_sumInew,
                                 double *gpu_sumRref, double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_energyRecipNew, int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  double sumRealOld = 0.0, sumImaginaryOld = 0.0;
  double sumRealNew = 0.0, sumImaginaryNew = 0.0;

  for (int p = 0; p < atomNumber; p++) {
    double dotProductOld =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_cx[p], gpu_cy[p], gpu_cz[p]);
    double dotProductNew =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_nx[p], gpu_ny[p], gpu_nz[p]);
    double oldsin, oldcos;
    sincos(dotProductOld, &oldsin, &oldcos);
    sumRealOld += (gpu_particleCharge[p] * oldcos);
    sumImaginaryOld += (gpu_particleCharge[p] * oldsin);
    double newsin, newcos;
    sincos(dotProductNew, &newsin, &newcos);
    sumRealNew += (gpu_particleCharge[p] * newcos);
    sumImaginaryNew += (gpu_particleCharge[p] * newsin);
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumRealOld + sumRealNew;
  gpu_sumInew[threadID] =
      gpu_sumIref[threadID] - sumImaginaryOld + sumImaginaryNew;

  gpu_energyRecipNew[threadID] =
      ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
        gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
       gpu_prefactRef[threadID]);
}

__global__ void ChangeLambdaMolReciprocalGPU(
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_kx, double *gpu_ky,
    double *gpu_kz, int atomNumber, double *gpu_particleCharge,
    double *gpu_sumRnew, double *gpu_sumInew, double *gpu_sumRref,
    double *gpu_sumIref, double *gpu_prefactRef, double *gpu_energyRecipNew,
    double lambdaCoef, int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  double sumRealNew = 0.0, sumImaginaryNew = 0.0;

  for (int p = 0; p < atomNumber; p++) {
    double dotProductNew =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_x[p], gpu_y[p], gpu_z[p]);
    double newsin, newcos;
    sincos(dotProductNew, &newsin, &newcos);
    sumRealNew += (gpu_particleCharge[p] * newcos);
    sumImaginaryNew += (gpu_particleCharge[p] * newsin);
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] + lambdaCoef * sumRealNew;
  gpu_sumInew[threadID] = gpu_sumIref[threadID] + lambdaCoef * sumImaginaryNew;

  gpu_energyRecipNew[threadID] =
      ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
        gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
       gpu_prefactRef[threadID]);
}

#endif
