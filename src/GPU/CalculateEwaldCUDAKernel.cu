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
  cudaMemset(vars->gpu_finalVal, 0, sizeof(double));
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  int threadsPerBlock = THREADS_PER_BLOCK_SM;
  int blocksPerGrid = std::min(static_cast<int>(imageSize), 5000);
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kx[box],
      vars->gpu_ky[box], vars->gpu_kz[box], vars->gpu_molCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], vars->gpu_prefact[box],
      vars->gpu_finalVal, atomNumber, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // // ReduceSum
  // DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
  // vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecip, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
  cudaMemset(vars->gpu_finalVal, 0, sizeof(double));
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  int threadsPerBlock = THREADS_PER_BLOCK_SM;
  int blocksPerGrid = std::min(static_cast<int>(imageSize), 5000);
  BoxReciprocalSumsGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_kxRef[box],
      vars->gpu_kyRef[box], vars->gpu_kzRef[box], vars->gpu_molCharge,
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], vars->gpu_prefactRef[box],
      vars->gpu_finalVal, atomNumber, imageSize);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  // DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
  // vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecip, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

// Note: This implementation assumes that this function is always called after
// BoxForce, so the coordinates have already been copied to the GPU. Otherwise,
// add cudaMemcpy calls to copy the coordinates to gpu_x, gpu_y, and gpu_z.
void CallBoxForceReciprocalGPU(
    VariablesCUDA *vars, XYZArray &atomForceRec, XYZArray &molForceRec,
    const std::vector<double> &particleCharge,
    const std::vector<int> &particleMol, const std::vector<int> &particleUsed,
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
      vars->gpu_startAtomIdx, vars->gpu_alpha, vars->gpu_alphaSq, constValue,
      imageSize, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], vars->gpu_x, vars->gpu_y, vars->gpu_z,
      vars->gpu_prefactRef[box], vars->gpu_sumRnew[box], vars->gpu_sumInew[box],
      vars->gpu_isFraction, vars->gpu_molIndex, vars->gpu_lambdaCoulomb,
      vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

void CallMolExchangeReciprocalGPU(VariablesCUDA *vars, uint imageSize, uint box,
                                  const std::vector<double> &molCharge,
                                  XYZArray const &molCoords,
                                  double &energyRecipNew, bool first_call) {
  // Calculate atom number -- exclude uncharged particles
  int atomNumber = molCharge.size();
  cudaMemcpy(vars->gpu_molCharge, &molCharge[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, &molCoords.x[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, &molCoords.y[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, &molCoords.z[0], atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = (imageSize + threadsPerBlock - 1) / threadsPerBlock;
  MolExchangeReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_kxRef[box], vars->gpu_kyRef[box], vars->gpu_kzRef[box],
      vars->gpu_prefactRef[box], vars->gpu_sumRref[box], vars->gpu_sumIref[box],
      vars->gpu_sumRnew[box], vars->gpu_sumInew[box], vars->gpu_molCharge,
      atomNumber, vars->gpu_x, vars->gpu_y, vars->gpu_z,
      vars->gpu_recipEnergies, imageSize, first_call);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_recipEnergies, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&energyRecipNew, vars->gpu_finalVal, sizeof(double),
             cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

// Calculate reciprocal term for lambdaNew and Old with same coordinates
void CallChangeLambdaMolReciprocalGPU(VariablesCUDA *vars,
                                      XYZArray const &coords,
                                      const std::vector<double> &molCharge,
                                      uint imageSize, double &energyRecipNew,
                                      const double lambdaCoef, uint box) {
  // Calculate atom number -- exclude uncharged particles
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

  int threadsPerBlock = THREADS_PER_BLOCK;
  int blocksPerGrid = (imageSize + threadsPerBlock - 1) / threadsPerBlock;
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
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif
}

__global__ void BoxReciprocalSumsGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_molCharge, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_prefact,
    double *__restrict__ gpu_finalVal, int atomNumber, int imageSize) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  // Specialize BlockReduce for a 1D block of threads of type double
  using BlockReduce = cub::BlockReduce<double, THREADS_PER_BLOCK_SM>;
  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage sumR_temp_storage;
  __shared__ typename BlockReduce::TempStorage sumI_temp_storage;

  for (int image = blockIdx.x; image < imageSize; image += gridDim.x) {
    double sumR = 0.0, sumI = 0.0;
#pragma unroll 8
    for (int particleID = threadIdx.x; particleID < atomNumber;
         particleID += THREADS_PER_BLOCK_SM) {
      double dot = DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image],
                                 gpu_x[particleID], gpu_y[particleID],
                                 gpu_z[particleID]);
      double dotsin, dotcos;
      sincos(dot, &dotsin, &dotcos);
      sumR += gpu_molCharge[particleID] * dotcos;
      sumI += gpu_molCharge[particleID] * dotsin;
    }
    __syncthreads();

    // Compute the block-wide sums for thread 0
    double aggregateR = BlockReduce(sumR_temp_storage).Sum(sumR);
    double aggregateI = BlockReduce(sumI_temp_storage).Sum(sumI);

    if (threadIdx.x == 0) {
      double recipEng = gpu_prefact[image];
      gpu_sumRnew[image] = aggregateR;
      gpu_sumInew[image] = aggregateI;
      recipEng *= aggregateR * aggregateR + aggregateI * aggregateI;
      atomicAdd(&gpu_finalVal[0], recipEng);
    }
  }
}

__global__ void BoxForceReciprocalGPU(
    double *__restrict__ gpu_aForceRecx, double *__restrict__ gpu_aForceRecy,
    double *__restrict__ gpu_aForceRecz, double *__restrict__ gpu_mForceRecx,
    double *__restrict__ gpu_mForceRecy, double *__restrict__ gpu_mForceRecz,
    const double *__restrict__ gpu_particleCharge,
    const int *__restrict__ gpu_particleMol,
    const int *__restrict__ gpu_particleUsed,
    const int *__restrict__ gpu_startAtomIdx, const double *gpu_alpha,
    const double *__restrict__ gpu_alphaSq, const double constValue,
    const int imageSize, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_prefact,
    const double *__restrict__ gpu_sumRnew,
    const double *__restrict__ gpu_sumInew,
    const bool *__restrict__ gpu_isFraction,
    const int *__restrict__ gpu_molIndex,
    const double *__restrict__ gpu_lambdaCoulomb,
    const double *__restrict__ gpu_cell_x,
    const double *__restrict__ gpu_cell_y,
    const double *__restrict__ gpu_cell_z,
    const double *__restrict__ gpu_Invcell_x,
    const double *__restrict__ gpu_Invcell_y,
    const double *__restrict__ gpu_Invcell_z,
    const int *__restrict__ gpu_nonOrth, double axx, double axy, double axz,
    int moveType, int box) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  __shared__ int shr_particleID, shr_moleculeID, shr_firstParticleID,
      shr_lastParticleID;
  __shared__ double shr_x, shr_y, shr_z, shr_lambdaCoefSq, shr_particleCharge,
      shr_fixed;

  if (threadIdx.x == 0) {
    // The particleID is the atom that corresponds to this particleUsed entry
    shr_particleID = gpu_particleUsed[blockIdx.x];
    shr_moleculeID = gpu_particleMol[shr_particleID];
    shr_x = gpu_x[shr_particleID];
    shr_y = gpu_y[shr_particleID];
    shr_z = gpu_z[shr_particleID];
    shr_firstParticleID = gpu_startAtomIdx[shr_moleculeID];
    shr_lastParticleID = gpu_startAtomIdx[shr_moleculeID + 1];
    double lambdaCoef = DeviceGetLambdaCoulomb(
        shr_moleculeID, box, gpu_isFraction, gpu_molIndex, gpu_lambdaCoulomb);
    shr_lambdaCoefSq = lambdaCoef * lambdaCoef;
    shr_particleCharge = gpu_particleCharge[shr_particleID];
    shr_fixed = 2.0 * lambdaCoef * shr_particleCharge;
  }
  __syncthreads();

  double forceX = 0.0, forceY = 0.0, forceZ = 0.0;

  // loop over images
#pragma unroll 8
  for (int image = threadIdx.x; image < imageSize; image += THREADS_PER_BLOCK) {
    double dotProduct = DotProductGPU(gpu_kx[image], gpu_ky[image],
                                      gpu_kz[image], shr_x, shr_y, shr_z);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    double factor = shr_fixed * gpu_prefact[image] *
                    (dotsin * gpu_sumRnew[image] - dotcos * gpu_sumInew[image]);
    forceX += factor * gpu_kx[image];
    forceY += factor * gpu_ky[image];
    forceZ += factor * gpu_kz[image];
  }

  // loop over other particles within the same molecule
  // pick the last warp -- the one most likely to exit the for loop early
  int lane = threadIdx.x - (THREADS_PER_BLOCK - 32);
  if (lane >= 0) {
    for (int otherParticle = shr_firstParticleID + lane;
         otherParticle < shr_lastParticleID; otherParticle += 32) {
      if (shr_particleID != otherParticle) {
        double distSq;
        double3 distVect;
        DeviceInRcut(distSq, distVect, gpu_x, gpu_y, gpu_z, shr_particleID,
                     otherParticle, axx, axy, axz, *gpu_nonOrth, gpu_cell_x,
                     gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                     gpu_Invcell_z);
        double qiqj =
            qqFactGPU * shr_particleCharge * gpu_particleCharge[otherParticle];
        double intraForce = qiqj * shr_lambdaCoefSq / distSq;
        double expVal = exp(-gpu_alphaSq[box] * distSq);
        double dist = sqrt(distSq);
        intraForce *= (erf(gpu_alpha[box] * dist) / dist) - constValue * expVal;
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
      gpu_aForceRecx[shr_particleID] = aggregateX;
      gpu_aForceRecy[shr_particleID] = aggregateY;
      gpu_aForceRecz[shr_particleID] = aggregateZ;
    } else if (moveType == mp::MPDISPLACE) {
      atomicAdd(&gpu_mForceRecx[shr_moleculeID], aggregateX);
      atomicAdd(&gpu_mForceRecy[shr_moleculeID], aggregateY);
      atomicAdd(&gpu_mForceRecz[shr_moleculeID], aggregateZ);
    }
  }
}

__global__ void SwapReciprocalGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    int atomNumber, const double *__restrict__ gpu_molCharge,
    double *__restrict__ gpu_sumRnew, double *__restrict__ gpu_sumInew,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, int imageSize) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  const int image = blockIdx.x * blockDim.x + threadIdx.x;
  if (image >= imageSize)
    return;

  double sumReal = gpu_sumRref[image], sumImag = gpu_sumIref[image];

#pragma unroll 6
  for (int p = 0; p < atomNumber; ++p) {
    double dotProduct =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image], gpu_x[p],
                      gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    // We negated the charge for molecule removals, so always add the charge
    sumImag += gpu_molCharge[p] * dotsin;
    sumReal += gpu_molCharge[p] * dotcos;
  }

  gpu_sumRnew[image] = sumReal;
  gpu_sumInew[image] = sumImag;

  gpu_recipEnergies[image] =
      (sumReal * sumReal + sumImag * sumImag) * gpu_prefactRef[image];
}

__global__ void __launch_bounds__(THREADS_PER_BLOCK) MolReciprocalGPU(
    const double *__restrict__ gpu_cx, const double *__restrict__ gpu_cy,
    const double *__restrict__ gpu_cz, const double *__restrict__ gpu_nx,
    const double *__restrict__ gpu_ny, const double *__restrict__ gpu_nz,
    const double *__restrict__ gpu_kx, const double *__restrict__ gpu_ky,
    const double *__restrict__ gpu_kz, int atomNumber,
    const double *__restrict__ gpu_molCharge, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, int imageSize) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  const int image = blockIdx.x * blockDim.x + threadIdx.x;
  if (image >= imageSize)
    return;

  double sumReal = gpu_sumRref[image], sumImag = gpu_sumIref[image];

#pragma unroll 4
  for (int p = 0; p < atomNumber; ++p) {
    double dotProductOld =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image], gpu_cx[p],
                      gpu_cy[p], gpu_cz[p]);
    double dotProductNew =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image], gpu_nx[p],
                      gpu_ny[p], gpu_nz[p]);
    double oldsin, oldcos;
    sincos(dotProductOld, &oldsin, &oldcos);
    double newsin, newcos;
    sincos(dotProductNew, &newsin, &newcos);
    sumImag += gpu_molCharge[p] * (newsin - oldsin);
    sumReal += gpu_molCharge[p] * (newcos - oldcos);
  }

  gpu_sumRnew[image] = sumReal;
  gpu_sumInew[image] = sumImag;

  gpu_recipEnergies[image] =
      (sumReal * sumReal + sumImag * sumImag) * gpu_prefactRef[image];
}

__global__ void MolExchangeReciprocalGPU(
    const double *__restrict__ gpu_kx, const double *__restrict__ gpu_ky,
    const double *__restrict__ gpu_kz,
    const double *__restrict__ gpu_prefactRef,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref, double *__restrict__ gpu_sumRnew,
    double *__restrict__ gpu_sumInew, const double *__restrict__ gpu_molCharge,
    int numChargedParticles, const double *__restrict__ gpu_x,
    const double *__restrict__ gpu_y, const double *__restrict__ gpu_z,
    double *__restrict__ gpu_recipEnergies, int imageSize, bool first_call) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  const int image = blockIdx.x * blockDim.x + threadIdx.x;
  if (image >= imageSize)
    return;

  double sumReal, sumImag;
  if (first_call) {
    sumReal = gpu_sumRref[image];
    sumImag = gpu_sumIref[image];
  } else {
    sumReal = gpu_sumRnew[image];
    sumImag = gpu_sumInew[image];
  }
#pragma unroll 6
  for (int p = 0; p < numChargedParticles; ++p) {
    double dotProduct =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image], gpu_x[p],
                      gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    sumImag += gpu_molCharge[p] * dotsin;
    sumReal += gpu_molCharge[p] * dotcos;
  }

  gpu_sumRnew[image] = sumReal;
  gpu_sumInew[image] = sumImag;

  gpu_recipEnergies[image] =
      (sumReal * sumReal + sumImag * sumImag) * gpu_prefactRef[image];
}

__global__ void ChangeLambdaMolReciprocalGPU(
    const double *__restrict__ gpu_x, const double *__restrict__ gpu_y,
    const double *__restrict__ gpu_z, const double *__restrict__ gpu_kx,
    const double *__restrict__ gpu_ky, const double *__restrict__ gpu_kz,
    int atomNumber, const double *__restrict__ gpu_molCharge,
    double *__restrict__ gpu_sumRnew, double *__restrict__ gpu_sumInew,
    const double *__restrict__ gpu_sumRref,
    const double *__restrict__ gpu_sumIref,
    const double *__restrict__ gpu_prefactRef,
    double *__restrict__ gpu_recipEnergies, double lambdaCoef, int imageSize) {
#if defined(NDEBUG) && (__CUDACC_VER_MAJOR__ >= 13)
  asm volatile(".pragma \"enable_smem_spilling\";");
#endif
  const int image = blockIdx.x * blockDim.x + threadIdx.x;
  if (image >= imageSize)
    return;

  double sumReal = gpu_sumRref[image], sumImag = gpu_sumIref[image];

#pragma unroll 6
  for (int p = 0; p < atomNumber; ++p) {
    double dotProduct =
        DotProductGPU(gpu_kx[image], gpu_ky[image], gpu_kz[image], gpu_x[p],
                      gpu_y[p], gpu_z[p]);
    double dotsin, dotcos;
    sincos(dotProduct, &dotsin, &dotcos);
    sumImag += lambdaCoef * gpu_molCharge[p] * dotsin;
    sumReal += lambdaCoef * gpu_molCharge[p] * dotcos;
  }

  gpu_sumRnew[image] = sumReal;
  gpu_sumInew[image] = sumImag;

  gpu_recipEnergies[image] =
      (sumReal * sumReal + sumImag * sumImag) * gpu_prefactRef[image];
}
#endif
