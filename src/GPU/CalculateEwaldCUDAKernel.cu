/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
                               real const *kx,
                               real const *ky,
                               real const *kz,
                               vector<real> particleCharge,
                               uint imageSize,
                               real *sumRnew,
                               real *sumInew,
                               real *prefact,
                               real *hsqr,
                               real &energyRecip,
                               uint box)
{
  real *gpu_particleCharge;
  real * gpu_energyRecip;
  real * gpu_final_energyRecip;
  int blocksPerGrid, threadsPerBlock;
  int atomNumber = coords.Count();

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_energyRecip, imageSize * sizeof(real));
  cudaMalloc((void**) &gpu_final_energyRecip, sizeof(real));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_kx[box], kx, imageSize * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ky[box], ky, imageSize * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_kz[box], kz, imageSize * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_prefact[box], prefact, imageSize * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_hsqr[box], hsqr, imageSize * sizeof(real),
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
             imageSize * sizeof(real),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box],
             imageSize * sizeof(real),
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
             sizeof(real), cudaMemcpyDeviceToHost);

  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecip);
  cudaFree(gpu_final_energyRecip);
  cudaFree(d_temp_storage);
}

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          vector<real> particleCharge,
                          uint imageSize,
                          real *sumRnew,
                          real *sumInew,
                          real &energyRecipNew,
                          uint box)
{
  // Calculate atom number
  int atomNumber = currentCoords.Count();
  int newCoordsNumber = newCoords.Count();
  real *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  real *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(real));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(real));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_nx, newCoords.x, newCoordsNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ny, newCoords.y, newCoordsNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nz, newCoords.z, newCoordsNumber * sizeof(real),
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
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(real),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(real),
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
             sizeof(real), cudaMemcpyDeviceToHost);


  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);
}

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           vector<real> particleCharge,
                           uint imageSize,
                           real *sumRnew,
                           real *sumInew,
                           int const insert,
                           real &energyRecipNew,
                           uint box)
{
  // Calculate atom number
  int atomNumber = coords.Count();
  // given coordinates
  real *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  real *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(real));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(real));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(real),
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

  // In the future maybe we could remove this for Nondebug?
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(real),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(real),
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
             sizeof(real), cudaMemcpyDeviceToHost);

  cudaFree(gpu_particleCharge);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);

}

__global__ void SwapReciprocalGPU(real *gpu_x, real *gpu_y, real *gpu_z,
                                  real *gpu_kx, real *gpu_ky, real *gpu_kz,
                                  int atomNumber,
                                  real *gpu_particleCharge,
                                  real *gpu_sumRnew,
                                  real *gpu_sumInew,
                                  real *gpu_sumRref,
                                  real *gpu_sumIref,
                                  real *gpu_prefactRef,
                                  int insert,
                                  real *gpu_energyRecipNew,
                                  int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  int p;
  real dotProduct = 0.0, sumReal = 0.0, sumImaginary = 0.0;

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

__global__ void MolReciprocalGPU(real *gpu_cx, real *gpu_cy, real *gpu_cz,
                                 real *gpu_nx, real *gpu_ny, real *gpu_nz,
                                 real *gpu_kx, real *gpu_ky, real *gpu_kz,
                                 int atomNumber,
                                 real *gpu_particleCharge,
                                 real *gpu_sumRnew,
                                 real *gpu_sumInew,
                                 real *gpu_sumRref,
                                 real *gpu_sumIref,
                                 real *gpu_prefactRef,
                                 real *gpu_energyRecipNew,
                                 int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int p;
  real dotProductOld = 0.0, dotProductNew = 0.0;
  real sumRealNew = 0.0, sumImaginaryNew = 0.0;
  real sumRealOld = 0.0, sumImaginaryOld = 0.0;

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

__global__ void BoxReciprocalSetupGPU(real *gpu_x,
                                      real *gpu_y,
                                      real *gpu_z,
                                      real *gpu_kx,
                                      real *gpu_ky,
                                      real *gpu_kz,
                                      real atomNumber,
                                      real *gpu_particleCharge,
                                      real *gpu_sumRnew,
                                      real *gpu_sumInew,
                                      int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int i;
  real dotP;

  gpu_sumRnew[threadID] = 0.0;
  gpu_sumInew[threadID] = 0.0;
  for(i = 0; i < atomNumber; i++) {
    dotP = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                         gpu_x[i], gpu_y[i], gpu_z[i]);
    gpu_sumRnew[threadID] += gpu_particleCharge[i] * cos(dotP);
    gpu_sumInew[threadID] += gpu_particleCharge[i] * sin(dotP);
  }
}

__global__ void BoxReciprocalGPU(real *gpu_prefact,
                                 real *gpu_sumRnew,
                                 real *gpu_sumInew,
                                 real *gpu_energyRecip,
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
