/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "CalculateEwaldCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "cub/cub.cuh"
#include <vector>

using namespace std;
using namespace cub;

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars,
                               XYZArray const &coords,
                               double const *kx,
                               double const *ky,
                               double const *kz,
                               vector<double> particleCharge,
                               uint imageSize,
                               double *sumRnew,
                               double *sumInew,
                               double *prefact,
                               double *hsqr,
                               double &energyRecip,
                               uint box)
{
  double *gpu_particleCharge;
  double * gpu_energyRecip;
  double * gpu_final_energyRecip;
  int blocksPerGrid, threadsPerBlock;
  int atomNumber = coords.Count();

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecip, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecip, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

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

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  BoxReciprocalSetupGPU <<< blocksPerGrid,
                        threadsPerBlock>>>(vars->gpu_x,
                            vars->gpu_y,
                            vars->gpu_z,
                            vars->gpu_kx[box],
                            vars->gpu_ky[box],
                            vars->gpu_kz[box],
                            atomNumber,
                            gpu_particleCharge,
                            vars->gpu_sumRnew[box],
                            vars->gpu_sumInew[box],
                            imageSize);

  BoxReciprocalGPU <<< blocksPerGrid, threadsPerBlock>>>(vars->gpu_prefact[box],
      vars->gpu_sumRnew[box],
      vars->gpu_sumInew[box],
      gpu_energyRecip,
      imageSize);

  // In the future maybe we could remove this for Nondebug?
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box],
             imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box],
             imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  cudaMemcpy(&energyRecip, gpu_final_energyRecip,
             sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecip);
  cudaFree(gpu_final_energyRecip);
  cudaFree(d_temp_storage);
}

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          vector<double> particleCharge,
                          uint imageSize,
                          double *sumRnew,
                          double *sumInew,
                          double &energyRecipNew,
                          uint box)
{
  // Calculate atom number
  int atomNumber = currentCoords.Count();
  int newCoordsNumber = newCoords.Count();
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

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
  MolReciprocalGPU <<< blocksPerGrid,
                   threadsPerBlock>>>(vars->gpu_x, vars->gpu_y, vars->gpu_z,
                                      vars->gpu_nx, vars->gpu_ny, vars->gpu_nz,
                                      vars->gpu_kxRef[box], vars->gpu_kyRef[box],
                                      vars->gpu_kzRef[box],
                                      atomNumber,
                                      gpu_particleCharge,
                                      vars->gpu_sumRnew[box],
                                      vars->gpu_sumInew[box],
                                      vars->gpu_sumRref[box],
                                      vars->gpu_sumIref[box],
                                      vars->gpu_prefactRef[box],
                                      gpu_energyRecipNew,
                                      imageSize);
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew,
             sizeof(double), cudaMemcpyDeviceToHost);


  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);
}

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           vector<double> particleCharge,
                           uint imageSize,
                           double *sumRnew,
                           double *sumInew,
                           int const insert,
                           double &energyRecipNew,
                           uint box)
{
  // Calculate atom number
  int atomNumber = coords.Count();
  // given coordinates
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  SwapReciprocalGPU <<< blocksPerGrid,
                    threadsPerBlock>>>(vars->gpu_x, vars->gpu_y, vars->gpu_z,
                                       vars->gpu_kxRef[box], vars->gpu_kyRef[box],
                                       vars->gpu_kzRef[box],
                                       atomNumber,
                                       gpu_particleCharge,
                                       vars->gpu_sumRnew[box],
                                       vars->gpu_sumInew[box],
                                       vars->gpu_sumRref[box],
                                       vars->gpu_sumIref[box],
                                       vars->gpu_prefactRef[box],
                                       insert,
                                       gpu_energyRecipNew,
                                       imageSize);

#ifndef NDEBUG
  // In the future maybe we could remove this for Nondebug?
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
#endif

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMalloc(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew,
             sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);

}

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                  double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                  int atomNumber,
                                  double *gpu_particleCharge,
                                  double *gpu_sumRnew,
                                  double *gpu_sumInew,
                                  double *gpu_sumRref,
                                  double *gpu_sumIref,
                                  double *gpu_prefactRef,
                                  int insert,
                                  double *gpu_energyRecipNew,
                                  int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  int p;
  double dotProduct = 0.0, sumReal = 0.0, sumImaginary = 0.0;

  for(p = 0; p < atomNumber; p++) {
    dotProduct = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                               gpu_kz[threadID], gpu_x[p], gpu_y[p], gpu_z[p]);
    sumReal += (gpu_particleCharge[p] * cos(dotProduct));
    sumImaginary += (gpu_particleCharge[p] * sin(dotProduct));
  }

  //If we insert the molecule to the box, we add the sum value.
  //Otherwise, we subtract the sum value
  if(insert) {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] + sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] + sumImaginary;
  } else {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] - sumImaginary;
  }

  gpu_energyRecipNew[threadID] = ((gpu_sumRnew[threadID] *
                                   gpu_sumRnew[threadID] +
                                   gpu_sumInew[threadID] *
                                   gpu_sumInew[threadID]) *
                                  gpu_prefactRef[threadID]);
}

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber,
                                 double *gpu_particleCharge,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_sumRref,
                                 double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_energyRecipNew,
                                 int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int p;
  double dotProductOld = 0.0, dotProductNew = 0.0;
  double sumRealNew = 0.0, sumImaginaryNew = 0.0;
  double sumRealOld = 0.0, sumImaginaryOld = 0.0;

  for(p = 0; p < atomNumber; p++) {
    dotProductOld = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                                  gpu_kz[threadID],
                                  gpu_cx[p], gpu_cy[p], gpu_cz[p]);
    dotProductNew = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                                  gpu_kz[threadID],
                                  gpu_nx[p], gpu_ny[p], gpu_nz[p]);
    sumRealNew += (gpu_particleCharge[p] * cos(dotProductNew));
    sumImaginaryNew += (gpu_particleCharge[p] * sin(dotProductNew));
    sumRealOld += (gpu_particleCharge[p] * cos(dotProductOld));
    sumImaginaryOld += (gpu_particleCharge[p] * sin(dotProductOld));
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumRealOld + sumRealNew;
  gpu_sumInew[threadID] = gpu_sumIref[threadID] - sumImaginaryOld +
                          sumImaginaryNew;

  gpu_energyRecipNew[threadID] = ((gpu_sumRnew[threadID] *
                                   gpu_sumRnew[threadID] +
                                   gpu_sumInew[threadID] *
                                   gpu_sumInew[threadID]) *
                                  gpu_prefactRef[threadID]);
}

__global__ void BoxReciprocalSetupGPU(double *gpu_x,
                                      double *gpu_y,
                                      double *gpu_z,
                                      double *gpu_kx,
                                      double *gpu_ky,
                                      double *gpu_kz,
                                      double atomNumber,
                                      double *gpu_particleCharge,
                                      double *gpu_sumRnew,
                                      double *gpu_sumInew,
                                      int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int i;
  double dotP;

  gpu_sumRnew[threadID] = 0.0;
  gpu_sumInew[threadID] = 0.0;
  for(i = 0; i < atomNumber; i++) {
    dotP = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                         gpu_x[i], gpu_y[i], gpu_z[i]);
    gpu_sumRnew[threadID] += gpu_particleCharge[i] * cos(dotP);
    gpu_sumInew[threadID] += gpu_particleCharge[i] * sin(dotP);
  }
}

__global__ void BoxReciprocalGPU(double *gpu_prefact,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_energyRecip,
                                 int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  gpu_energyRecip[threadID] = ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
                                gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
                               gpu_prefact[threadID]);
}

#endif
