#include "CalculateEwaldCUDAKernel.h"

#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDA.h"
#include "CalculateMinImageCUDA.h"

void CallBoxReciprocalSetupGPU(XYZArray const & coords,
			       double const *kx,
			       double const *ky,
			       double const *kz,
			       vector<double> particleCharge,
			       uint imageSize,
			       double *sumRnew,
			       double *sumInew)
{
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_kx, *gpu_ky, *gpu_kz;
  double *gpu_particleCharge;
  double *gpu_sumRnew, *gpu_sumInew;
  double *gpu_imageSize;
  int start, length;
  int i = 0;
  int blocksPerGrid, threadsPerBlock;
  int atomNumber = coords.Count();

  cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size * sizeof(double));
  cudaMalloc((void**) &gpu_sumRnew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumInew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kx, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_ky, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kz, imageSize * sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0], 
	     particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_x, coords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, coords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, coords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kx, kx, imageSize * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_ky, ky, imageSize * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kz, kz, imageSize * sizeof(double), cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize/threadsPerBlock) + 1;

  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_sumRnew);
  cudaFree(gpu_sumInew);
  cudaFree(gpu_kx);
  cudaFree(gpu_ky);
  cudaFree(gpu_kz);
}

void CallBoxReciprocalGPU(double * prefact,
			  double * sumRnew,
			  double * sumInew,
			  int imageSize)
{
  double * gpu_sumRnew, *gpu_sumInew;
  double * gpu_prefact;
  double * gpu_energyRecip;
  int blockPerGrid, threadsPerBlock;
  
  cudaMalloc((void**) &gpu_sumRnew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumInew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_prefact, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecip, imageSize * sizeof(double));

  cudaMemcpy(gpu_prefact, prefact, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumRnew, sumRnew, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumInew, sumInew, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  
  threadsPerBlock = 256;
  blocksPerGrid = (imageSize/threadsPerBlock) + 1;
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_prefact,
						       gpu_sumRnew,
						       gpu_sumInew,
						       gpu_energyRecip,
						       imageSize);

  cudaFree(gpu_sumRnew);
  cudaFree(gpu_sumInew);
  cudaFree(gpu_prefact);
}

__global__ void BoxReciprocalSetupGPU(double * gpu_x,
				      double * gpu_y,
				      dobule * gpu_z,
				      double * gpu_kx,
				      double * gpu_ky,
				      double * gpu_kz,
				      double atomNumber,
				      double * gpu_particleCharge,
				      double * gpu_sumRnew,
				      double * gpu_sumInew,
				      double imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID>=imageSize)
    return;
  int i;
  double dotP;
  
  sumRnew[threadiD] = 0.0;
  sumInew[threadID] = 0.0;
  for(i=0; i<atomNumber; i++)
  {
    dotP = DotProduct(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
		      gpu_x[threadID], gpu_y[threadID], gpu_z[threadID]);
    gpu_sumRnew[threadID] += gpu_particleCharge[threadID] * cos(dotP);
    gpu_sumInew[threadID] += gpu_particleCharge[threadID] * sin(dotP);
  }

  // reduction
}

__global__ void BoxReciprocalGPU(double *gpu_prefact,
				 double *gpu_sumRnew,
				 double *gpu_sumInew,
				 double *gpu_energyRecip,
				 int mageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID>=imageSize)
    return;
  
  gpu_energyRecip[threadID] = gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
    gpu_sumInew[threadID] * gpu_sumInew[threadID] * gpu_prefact[threadID];
}

__device__ double DotProduct(double kx, double ky, double kx, 
			     double x, double y, double z)
{
  return (kx * x + ky * y + kz * z);
}

#endif