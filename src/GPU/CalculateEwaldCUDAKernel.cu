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
#include "MoveSettings.h"
#include "cub/cub.cuh"

using namespace cub;

const int THREADS_PER_BLOCK = 128;
const int THREADS_PER_BLOCK_SM = 64;

// Use this function when calculating the reciprocal terms
// for a new volume. Such as a change in box dimensions.
void CallBoxReciprocalSetupGPU(VariablesCUDA *vars, XYZArray const &coords,
                               double const *kx, double const *ky,
                               double const *kz,
                               const std::vector<double> &molCharge,
                               uint imageSize, double *prefact, double *hsqr,
                               double &energyRecip, uint box) {
  int atomNumber = molCharge.size();

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_kx[box], kx, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ky[box], ky, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_kz[box], kz, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_prefact[box], prefact, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_hsqr[box], hsqr, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  int threadsPerBlock = THREADS_PER_BLOCK_SM;
  int blocksPerGrid = imageSize;
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kx[box],
      vars->gpu_ky[box], vars->gpu_kz[box], atomNumber, vars->gpu_molCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box]);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // Fewer blocks are needed since we are doing one computation per image
  threadsPerBlock = THREADS_PER_BLOCK;
  blocksPerGrid = (imageSize + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_prefact[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_recipEnergies, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecip, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

// Use this function when calculating the reciprocal terms
// with the current volume.
void CallBoxReciprocalSumsGPU(VariablesCUDA *vars, XYZArray const &coords,
                              const std::vector<double> &molCharge,
                              uint imageSize, double &energyRecip, uint box) {
  int atomNumber = molCharge.size();

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  int threadsPerBlock = THREADS_PER_BLOCK_SM;
  int blocksPerGrid = imageSize;
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], atomNumber,
      vars->gpu_molCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box]);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // Fewer blocks are needed since we are doing one computation per image
  threadsPerBlock = THREADS_PER_BLOCK;
  blocksPerGrid = (imageSize + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_prefactRef[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_recipEnergies, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecip, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

__global__ void BoxReciprocalSumsGPU(double *gpu_x, double *gpu_y,
                                     double *gpu_z, double *gpu_kx,
                                     double *gpu_ky, double *gpu_kz,
                                     int atomNumber, double *gpu_molCharge,
                                     double *gpu_sumRnew, double *gpu_sumInew) {
  int image = blockIdx.x;
  double sumR = 0.0, sumI = 0.0;
#pragma unroll 8
  for (int particleID = threadIdx.x; particleID < atomNumber;
       particleID += THREADS_PER_BLOCK_SM) {
    double dot =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image],
                      gpu_x[particleID], gpu_y[particleID], gpu_z[particleID]);
    double dotsin, dotcos;
    sincos(dot, &dotsin, &dotcos);
    sumR += gpu_molCharge[particleID] * dotcos;
    sumI += gpu_molCharge[particleID] * dotsin;
  }
  __syncthreads();

  // Specialize BlockReduce for a 1D block of threads of type double
  using BlockReduce = cub::BlockReduce<double, THREADS_PER_BLOCK_SM>;

  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage sumR_temp_storage;
  __shared__ typename BlockReduce::TempStorage sumI_temp_storage;

  // Compute the block-wide sums for thread 0
  double aggregateR = BlockReduce(sumR_temp_storage).Sum(sumR);
  double aggregateI = BlockReduce(sumI_temp_storage).Sum(sumI);

  if (threadIdx.x == 0) {
    gpu_sumRnew[image] = aggregateR;
    gpu_sumInew[image] = aggregateI;
  }
}

__global__ void BoxReciprocalGPU(double *gpu_prefact, double *gpu_sumRnew,
                                 double *gpu_sumInew, double *gpu_recipEnergies,
                                 int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  gpu_recipEnergies[threadID] =
      (gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
       gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
      gpu_prefact[threadID];
}

void CallMolReciprocalGPU(VariablesCUDA *vars, XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          const std::vector<double> &molCharge, uint imageSize,
                          double &energyRecipNew, uint box) {
  // Calculate number of coordinates -- exclude uncharged particles
  int CoordsNumber = molCharge.size();

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, currentCoords.x, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_nx, newCoords.x, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ny, newCoords.y, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nz, newCoords.z, CoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = (imageSize + threadsPerBlock - 1) / threadsPerBlock;
  MolReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_nx, vars->gpu_ny,
      vars->gpu_nz, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], CoordsNumber, vars->gpu_molCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], vars->gpu_sumRref[box],
      vars->gpu_sumIref[box], vars->gpu_prefactRef[box],
      vars->gpu_recipEnergies, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecipNew, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

// Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &molCharge,
                                      uint imageSize, double &energyRecipNew,
                                      const double lambdaCoef, uint box) {
  // Calculate atom number -- exclude uncharged particles
  int atomNumber = molCharge.size();
  int blocksPerGrid, threadsPerBlock;

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = THREADS_PER_BLOCK;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  ChangeLambdaMolReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], atomNumber,
      vars->gpu_molCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_sumRref[box], vars->gpu_sumIref[box], vars->gpu_prefactRef[box],
      vars->gpu_recipEnergies, lambdaCoef, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecipNew, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

void CallSwapReciprocalGPU(VariablesCUDA *vars, XYZArray const &coords,
                           const std::vector<double> &molCharge, uint imageSize,
                           double &energyRecipNew, uint box) {
  // Calculate atom number -- exclude uncharged particles
  int atomNumber = molCharge.size();

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_nx, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ny, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nz, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = (imageSize + threadsPerBlock - 1) / threadsPerBlock;
  SwapReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_nx, vars->gpu_ny, vars->gpu_nz, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], molCharge.size(),
      vars->gpu_molCharge, vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_sumRref[box], vars->gpu_sumIref[box], vars->gpu_prefactRef[box],
      vars->gpu_recipEnergies, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecipNew, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars, uint imageSize, uint box,
                                  const std::vector<double> &molCharge,
                                  double &energyRecipNew,
                                  XYZArray const &molCoords) {
  // Calculate atom number -- exclude uncharged particles
  int atomNumber = molCharge.size();

  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, &molCoords.x[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, &molCoords.y[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, &molCoords.z[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = (imageSize + threadsPerBlock - 1) / threadsPerBlock;
  MolExchangeReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      imageSize, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], vars->gpu_prefactRef[box], vars->gpu_sumRnew[box],
      vars->gpu_sumInew[box], vars->gpu_molCharge, atomNumber, vars->gpu_x,
      vars->gpu_y, vars->gpu_z, vars->gpu_recipEnergies);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecipNew, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
}

// Note: This implementation assumes that this function is always called after
// BoxForce, so the coordinates have already been copied to the GPU. Otherwise,
// add cudaMemcpy calls to copy the coordinates to gpu_x, gpu_y, and gpu_z.
void CallBoxForceReciprocalGPU(
    VariablesCUDA *vars, XYZArray &atomForceRec, XYZArray &molForceRec,
    const std::vector<double> &particleCharge,
    const std::vector<int> &particleMol, const std::vector<int> &particleUsed,
    const std::vector<int> &startMol, const std::vector<int> &lengthMol,
    double constValue, uint imageSize, XYZArray const &molCoords,
    BoxDimensions const &boxAxes, int moveType, int box) {
  int atomCount = atomForceRec.Count();
  int molCount = molForceRec.Count();

  // calculate block and grid sizes
  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = particleUsed.size();

  if (moveType == mp::MPDISPLACE) {
    cudaMemset(vars->gpu_mForceRecx, 0, molCount * sizeof(double));
    cudaMemset(vars->gpu_mForceRecy, 0, molCount * sizeof(double));
    cudaMemset(vars->gpu_mForceRecz, 0, molCount * sizeof(double));
  }
  cudaMemcpy(vars->gpu_particleUsed, &particleUsed[0],
             sizeof(int) * particleUsed.size(), cudaMemcpyHostToDevice);

#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  BoxForceReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_aForceRecx, vars->gpu_aForceRecy, vars->gpu_aForceRecz,
      vars->gpu_mForceRecx, vars->gpu_mForceRecy, vars->gpu_mForceRecz,
      vars->gpu_particleCharge, vars->gpu_particleMol, vars->gpu_particleUsed,
      vars->gpu_startMol, vars->gpu_lengthMol, vars->gpu_alpha,
      vars->gpu_alphaSq, constValue, imageSize, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], vars->gpu_x, vars->gpu_y,
      vars->gpu_z, vars->gpu_prefactRef[box], vars->gpu_sumRnew[box],
      vars->gpu_sumInew[box], vars->gpu_isFraction, vars->gpu_molIndex,
      vars->gpu_lambdaCoulomb, vars->gpu_cell_x[box], vars->gpu_cell_y[box],
      vars->gpu_cell_z[box], vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box], vars->gpu_nonOrth, boxAxes.GetAxis(box).x,
      boxAxes.GetAxis(box).y, boxAxes.GetAxis(box).z, moveType, box);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  if (moveType == mp::MPROTATE) {
    cudaMemcpy(atomForceRec.x, vars->gpu_aForceRecx, sizeof(double) * atomCount,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(atomForceRec.y, vars->gpu_aForceRecy, sizeof(double) * atomCount,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(atomForceRec.z, vars->gpu_aForceRecz, sizeof(double) * atomCount,
               cudaMemcpyDeviceToHost);
  } else if (moveType == mp::MPDISPLACE) {
    cudaMemcpy(molForceRec.x, vars->gpu_mForceRecx, sizeof(double) * molCount,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(molForceRec.y, vars->gpu_mForceRecy, sizeof(double) * molCount,
               cudaMemcpyDeviceToHost);
    cudaMemcpy(molForceRec.z, vars->gpu_mForceRecz, sizeof(double) * molCount,
               cudaMemcpyDeviceToHost);
  }
}

__global__ void BoxForceReciprocalGPU(
    double *gpu_aForceRecx, double *gpu_aForceRecy, double *gpu_aForceRecz,
    double *gpu_mForceRecx, double *gpu_mForceRecy, double *gpu_mForceRecz,
    double *gpu_particleCharge, int *gpu_particleMol,
    const int *gpu_particleUsed, int *gpu_startMol, int *gpu_lengthMol,
    double *gpu_alpha, double *gpu_alphaSq, double constValue, int imageSize,
    double *gpu_kx, double *gpu_ky, double *gpu_kz, double *gpu_x,
    double *gpu_y, double *gpu_z, double *gpu_prefact, double *gpu_sumRnew,
    double *gpu_sumInew, bool *gpu_isFraction, int *gpu_molIndex,
    double *gpu_lambdaCoulomb, double *gpu_cell_x, double *gpu_cell_y,
    double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
    double *gpu_Invcell_z, int *gpu_nonOrth, double axx, double axy, double axz,
    int moveType, int box) {

  __shared__ int particleID, moleculeID;
  __shared__ double x, y, z, lambdaCoef, fixed;

  if (threadIdx.x == 0) {
    // The particleID is the atom that corresponds to this particleUsed entry
    particleID = gpu_particleUsed[blockIdx.x];
    moleculeID = gpu_particleMol[particleID];
    x = gpu_x[particleID];
    y = gpu_y[particleID];
    z = gpu_z[particleID];
    lambdaCoef = DeviceGetLambdaCoulomb(moleculeID, box, gpu_isFraction,
                                        gpu_molIndex, gpu_lambdaCoulomb);
    fixed = 2.0 * lambdaCoef * gpu_particleCharge[particleID];
  }
  __syncthreads();

  double forceX = 0.0, forceY = 0.0, forceZ = 0.0;

  // loop over images
  for (int image = threadIdx.x; image < imageSize; image += THREADS_PER_BLOCK) {
    double dot = x * gpu_kx[image] + y * gpu_ky[image] + z * gpu_kz[image];
    double dotsin, dotcos;
    sincos(dot, &dotsin, &dotcos);
    double factor = fixed * gpu_prefact[image] *
                    (dotsin * gpu_sumRnew[image] - dotcos * gpu_sumInew[image]);
    forceX += factor * gpu_kx[image];
    forceY += factor * gpu_ky[image];
    forceZ += factor * gpu_kz[image];
  }

  // loop over other particles within the same molecule
  // Pick the thread most likely to exit the for loop early
  if (threadIdx.x == THREADS_PER_BLOCK - 1) {
    double intraForce = 0.0, distSq = 0.0, dist = 0.0;
    double3 distVect;
    int lastParticleWithinSameMolecule =
        gpu_startMol[particleID] + gpu_lengthMol[particleID];
    for (int otherParticle = gpu_startMol[particleID];
         otherParticle < lastParticleWithinSameMolecule; ++otherParticle) {
      if (particleID != otherParticle) {
        DeviceInRcut(distSq, distVect, gpu_x, gpu_y, gpu_z, particleID,
                     otherParticle, axx, axy, axz, *gpu_nonOrth, gpu_cell_x,
                     gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                     gpu_Invcell_z);
        dist = sqrt(distSq);

        double expConstValue = exp(-1.0 * gpu_alphaSq[box] * distSq);
        double qiqj = gpu_particleCharge[particleID] *
                      gpu_particleCharge[otherParticle] * qqFactGPU;
        intraForce = qiqj * lambdaCoef * lambdaCoef / distSq;
        intraForce *=
            (erf(gpu_alpha[box] * dist) / dist) - constValue * expConstValue;
        forceX -= intraForce * distVect.x;
        forceY -= intraForce * distVect.y;
        forceZ -= intraForce * distVect.z;
      }
    }
  }
  __syncthreads();

  // Specialize BlockReduce for a 1D block of threads of type double
  using BlockReduce = cub::BlockReduce<double, THREADS_PER_BLOCK>;

  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage forceX_temp_storage;
  __shared__ typename BlockReduce::TempStorage forceY_temp_storage;
  __shared__ typename BlockReduce::TempStorage forceZ_temp_storage;

  // Compute the block-wide sum for thread 0
  double aggregateX = BlockReduce(forceX_temp_storage).Sum(forceX);
  double aggregateY = BlockReduce(forceY_temp_storage).Sum(forceY);
  double aggregateZ = BlockReduce(forceZ_temp_storage).Sum(forceZ);

  if (threadIdx.x == 0) {
    if (moveType == mp::MPROTATE) {
      gpu_aForceRecx[particleID] = aggregateX;
      gpu_aForceRecy[particleID] = aggregateY;
      gpu_aForceRecz[particleID] = aggregateZ;
    } else if (moveType == mp::MPDISPLACE) {
      atomicAdd(&gpu_mForceRecx[moleculeID], aggregateX);
      atomicAdd(&gpu_mForceRecy[moleculeID], aggregateY);
      atomicAdd(&gpu_mForceRecz[moleculeID], aggregateZ);
    }
  }
}

__global__ void SwapReciprocalGPU(
    const double *gpu_x, const double *gpu_y, const double *gpu_z,
    const double *gpu_kx, const double *gpu_ky, const double *gpu_kz,
    int atomNumber, const double *gpu_molCharge, double *gpu_sumRnew,
    double *gpu_sumInew, const double *gpu_sumRref, const double *gpu_sumIref,
    const double *gpu_prefactRef, double *gpu_recipEnergies, int imageSize) {

  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  double sumReal = gpu_sumRref[threadID], sumImaginary = gpu_sumIref[threadID];

#pragma unroll 4
  for (int p = 0; p < atomNumber; ++p) {
    double dotProduct =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_x[p], gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    // We negated the charge for molecule removals, so always add the charge
    sumReal += gpu_molCharge[p] * dotcos;
    sumImaginary += gpu_molCharge[p] * dotsin;
  }

  gpu_sumRnew[threadID] = sumReal;
  gpu_sumInew[threadID] = sumImaginary;

  gpu_recipEnergies[threadID] =
      ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
        gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
       gpu_prefactRef[threadID]);
}

__global__ void MolExchangeReciprocalGPU(
    int imageSize, double *gpu_kx, double *gpu_ky, double *gpu_kz,
    double *gpu_prefactRef, double *gpu_sumRnew, double *gpu_sumInew,
    double *gpu_molCharge, int numChargedParticles, double *gpu_x,
    double *gpu_y, double *gpu_z, double *gpu_recipEnergies) {
  int imageID = blockIdx.x * blockDim.x + threadIdx.x;
  if (imageID >= imageSize)
    return;

  double sumRnew = gpu_sumRnew[imageID], sumInew = gpu_sumInew[imageID];
#pragma unroll 6
  for (int p = 0; p < numChargedParticles; ++p) {
    double dotProduct =
        DotProductGPU(gpu_kx[imageID], gpu_ky[imageID], gpu_kz[imageID],
                      gpu_x[p], gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    sumRnew += gpu_molCharge[p] * dotcos;
    sumInew += gpu_molCharge[p] * dotsin;
  }

  gpu_sumRnew[imageID] = sumRnew;
  gpu_sumInew[imageID] = sumInew;
  gpu_recipEnergies[imageID] =
      (sumRnew * sumRnew + sumInew * sumInew) * gpu_prefactRef[imageID];
}

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber, double *gpu_molCharge,
                                 double *gpu_sumRnew, double *gpu_sumInew,
                                 double *gpu_sumRref, double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_recipEnergies, int imageSize) {
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if (threadID >= imageSize)
    return;

  double sumReal = 0.0, sumImaginary = 0.0;

#pragma unroll 4
  for (int p = 0; p < atomNumber; ++p) {
    double dotProductOld =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_cx[p], gpu_cy[p], gpu_cz[p]);
    double dotProductNew =
        DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                      gpu_nx[p], gpu_ny[p], gpu_nz[p]);
    double oldsin, oldcos;
    sincos(dotProductOld, &oldsin, &oldcos);
    sumReal -= gpu_molCharge[p] * oldcos;
    sumImaginary -= gpu_molCharge[p] * oldsin;
    double newsin, newcos;
    sincos(dotProductNew, &newsin, &newcos);
    sumReal += gpu_molCharge[p] * newcos;
    sumImaginary += gpu_molCharge[p] * newsin;
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] + sumReal;
  gpu_sumInew[threadID] = gpu_sumIref[threadID] + sumImaginary;

  gpu_recipEnergies[threadID] =
      ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
        gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
       gpu_prefactRef[threadID]);
}

__global__ void ChangeLambdaMolReciprocalGPU(
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_kx, double *gpu_ky,
    double *gpu_kz, int atomNumber, double *gpu_molCharge, double *gpu_sumRnew,
    double *gpu_sumInew, double *gpu_sumRref, double *gpu_sumIref,
    double *gpu_prefactRef, double *gpu_recipEnergies, double lambdaCoef,
    int imageSize) {
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
    sumRealNew += gpu_molCharge[p] * newcos;
    sumImaginaryNew += gpu_molCharge[p] * newsin;
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] + lambdaCoef * sumRealNew;
  gpu_sumInew[threadID] = gpu_sumIref[threadID] + lambdaCoef * sumImaginaryNew;

  gpu_recipEnergies[threadID] =
      (gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
       gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
      gpu_prefactRef[threadID];
}
#endif
