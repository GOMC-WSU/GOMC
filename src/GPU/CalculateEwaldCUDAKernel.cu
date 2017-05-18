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

void CallBoxReciprocalSetupGPU(XYZArray const &coords,
			       double const *kx,
			       double const *ky,
			       double const *kz,
			       vector<double> particleCharge,
			       uint imageSize,
			       double *sumRnew,
			       double *sumInew,
			       double *prefact,
			       double &energyRecip)
{
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_kx, *gpu_ky, *gpu_kz;
  double *gpu_particleCharge;
  double *gpu_sumRnew, *gpu_sumInew;
  double * gpu_prefact;
  double * gpu_energyRecip;
  double * gpu_final_energyRecip;
  int blocksPerGrid, threadsPerBlock;
  int atomNumber = coords.Count();

  cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_sumRnew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumInew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kx, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_ky, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kz, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_prefact, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecip, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecip, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
	     particleCharge.size() * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_x, coords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, coords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, coords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kx, kx, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_ky, ky, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kz, kz, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_prefact, prefact, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize/threadsPerBlock) + 1;
  BoxReciprocalSetupGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_x, gpu_y, 
							    gpu_z, gpu_kx, 
							    gpu_ky, gpu_kz, 
							    atomNumber, 
							    gpu_particleCharge,
							    gpu_sumRnew,
							    gpu_sumInew, 
							    imageSize);
  BoxReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_prefact,
						       gpu_sumRnew,
						       gpu_sumInew,
						       gpu_energyRecip,
						       imageSize);

  cudaMemcpy(sumRnew, gpu_sumRnew, imageSize * sizeof(double), 
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, gpu_sumInew, imageSize * sizeof(double),
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

  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_sumRnew);
  cudaFree(gpu_sumInew);
  cudaFree(gpu_kx);
  cudaFree(gpu_ky);
  cudaFree(gpu_kz);
  cudaFree(gpu_prefact);
  cudaFree(gpu_energyRecip);
  cudaFree(gpu_final_energyRecip);
  cudaFree(d_temp_storage);
}

void CallMolReciprocalGPU(XYZArray const &currentCoords,
			  XYZArray const &newCoords,
			  double const *kx,
			  double const *ky,
			  double const *kz,
			  vector<double> particleCharge,
			  uint imageSize,
			  double const *sumRref,
			  double const *sumIref,
			  double *sumRnew,
			  double *sumInew,
			  double const *prefactRef,
			  double &energyRecipNew)
{
  // Calculate atom number
  int atomNumber = currentCoords.Count();
  // Current coordinates
  double *gpu_cx, *gpu_cy, *gpu_cz;
  // New coordinates
  double *gpu_nx, *gpu_ny, *gpu_nz;
  double *gpu_kx, *gpu_ky, *gpu_kz;
  double *gpu_particleCharge;
  double *gpu_sumRnew, *gpu_sumInew;
  double *gpu_sumRref, *gpu_sumIref;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_prefactRef;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_cx, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_cy, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_cz, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_nx, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_ny, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_nz, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_sumRnew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumInew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumRref, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumIref, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kx, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_ky, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kz, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_prefactRef, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
	     particleCharge.size() * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_cx, currentCoords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cy, currentCoords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cz, currentCoords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_nx, newCoords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_ny, newCoords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_nz, newCoords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_kx, kx, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_ky, ky, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kz, kz, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_prefactRef, prefactRef, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumRref, sumRref, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumIref, sumIref, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize/threadsPerBlock) + 1;
  MolReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_cx, gpu_cy, gpu_cz,
						       gpu_nx, gpu_ny, gpu_nz,
						       gpu_kx, gpu_ky, gpu_kz, 
						       atomNumber, 
						       gpu_particleCharge,
						       gpu_sumRnew,
						       gpu_sumInew,
						       gpu_sumRref,
						       gpu_sumIref,
						       gpu_prefactRef,
						       gpu_energyRecipNew,
						       imageSize);
  cudaMemcpy(sumRnew, gpu_sumRnew, imageSize * sizeof(double),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, gpu_sumInew, imageSize * sizeof(double),
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

  
  cudaFree(gpu_cx);
  cudaFree(gpu_cy);
  cudaFree(gpu_cz);
  cudaFree(gpu_nx);
  cudaFree(gpu_ny);
  cudaFree(gpu_nz);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_sumRref);
  cudaFree(gpu_sumIref);
  cudaFree(gpu_sumRnew);
  cudaFree(gpu_sumInew);
  cudaFree(gpu_kx);
  cudaFree(gpu_ky);
  cudaFree(gpu_kz);
  cudaFree(gpu_prefactRef);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);
}

void CallSwapReciprocalGPU(XYZArray const &coords,
			   double const *kx,
			   double const *ky,
			   double const *kz,
			   vector<double> particleCharge,
			   uint imageSize,
			   double const *sumRref,
			   double const *sumIref,
			   double *sumRnew,
			   double *sumInew,
			   double const *prefactRef,
			   int const insert,
			   double &energyRecipNew)
{
  // Calculate atom number
  int atomNumber = coords.Count();
  // given coordinates
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_kx, *gpu_ky, *gpu_kz;
  double *gpu_particleCharge;
  double *gpu_sumRnew, *gpu_sumInew;
  double *gpu_sumRref, *gpu_sumIref;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_prefactRef;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_sumRnew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumInew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumRref, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_sumIref, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kx, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_ky, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_kz, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_prefactRef, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  cudaMalloc((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
	     particleCharge.size() * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_x, coords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, coords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, coords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_kx, kx, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_ky, ky, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_kz, kz, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_prefactRef, prefactRef, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumRref, sumRref, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_sumIref, sumIref, imageSize * sizeof(double),
	     cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize/threadsPerBlock) + 1;
  SwapReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_x, gpu_y, gpu_z,
							gpu_kx, gpu_ky, gpu_kz, 
							atomNumber, 
							gpu_particleCharge,
							gpu_sumRnew,
							gpu_sumInew,
							gpu_sumRref,
							gpu_sumIref,
							gpu_prefactRef,
							insert,
							gpu_energyRecipNew,
							imageSize);

  cudaMemcpy(sumRnew, gpu_sumRnew, imageSize * sizeof(double),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, gpu_sumInew, imageSize * sizeof(double),
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

  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_sumRref);
  cudaFree(gpu_sumIref);
  cudaFree(gpu_sumRnew);
  cudaFree(gpu_sumInew);
  cudaFree(gpu_kx);
  cudaFree(gpu_ky);
  cudaFree(gpu_kz);
  cudaFree(gpu_prefactRef);
  cudaFree(gpu_energyRecipNew);
  cudaFree(gpu_final_energyRecipNew);
  cudaFree(d_temp_storage);

}

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
				  double *gpu_kx, double *gpu_ky,double *gpu_kz,
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
  if(threadID>=imageSize)
    return;

  int p;
  double dotProduct = 0.0, sumReal = 0.0, sumImaginary = 0.0;

  for(p = 0; p < atomNumber; p++)
  {
    dotProduct = DotProduct(gpu_kx[threadID], gpu_ky[threadID],gpu_kz[threadID],
			    gpu_x[p], gpu_y[p], gpu_z[p]);
    sumReal += (gpu_particleCharge[p] * cos(dotProduct));
    sumImaginary += (gpu_particleCharge[p] * sin(dotProduct));
  }

  //If we insert the molecule to the box, we add the sum value.
  //Otherwise, we subtract the sum value
  if(insert)
  {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] + sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] + sumImaginary;
  }
  else
  {
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
  if(threadID>=imageSize)
    return;
  int p;
  double dotProductOld = 0.0, dotProductNew = 0.0;
  double sumRealNew = 0.0, sumImaginaryNew = 0.0;
  double sumRealOld = 0.0, sumImaginaryOld = 0.0; 
  
  for(p = 0; p < atomNumber; p++)
  {
    dotProductOld = DotProduct(gpu_kx[threadID], gpu_ky[threadID], 
			       gpu_kz[threadID],
			       gpu_cx[p], gpu_cy[p], gpu_cz[p]);
    dotProductNew = DotProduct(gpu_kx[threadID], gpu_ky[threadID], 
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
  if(threadID>=imageSize)
    return;
  int i;
  double dotP;
  
  gpu_sumRnew[threadID] = 0.0;
  gpu_sumInew[threadID] = 0.0;
  for(i = 0; i<atomNumber; i++)
  {
    dotP = DotProduct(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
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
  if(threadID>=imageSize)
    return;
  
  gpu_energyRecip[threadID] = ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
				gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
			       gpu_prefact[threadID]);
}

__device__ double DotProduct(double kx, double ky, double kz, 
			     double x, double y, double z)
{
  return (kx * x + ky * y + kz * z);
}

#endif
