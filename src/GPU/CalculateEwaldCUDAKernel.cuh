#ifndef CALCULATE_EWALD_CUDA_KERNEL
#define CALCULATE_EWALD_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "XYZArray.h"

using namespace std;

void CallBoxReciprocalSetupGPU(XYZArray const & coords,
			       double const *kx,
			       double const *ky,
			       double const *kz,
			       vector<double> particleCharge,
			       uint imageSize,
			       double *sumRnew,
			       double *sumInew);

__global__ void BoxReciprocalSetupGPU(double * gpu_x,
				      double * gpu_y,
				      double * gpu_z,
				      double * gpu_kx,
				      double * gpu_ky,
				      double * gpu_kz,
				      double atomNumber,
				      double * gpu_particleCharge,
				      double * gpu_sumRnew,
				      double * gpu_sumInew,
				      double imageSize);

void CallBoxReciprocalGPU(double * prefact,
			  double * sumRnew,
			  double * sumInew,
			  double & energyRecip,
			  int imageSize);

__global__ void BoxReciprocalGPU(double *gpu_prefact,
				 double *gpu_sumRnew,
				 double *gpu_sumInew,
				 double *gpu_energyRecip,
				 int mageSize);

__device__ double DotProduct(double kx, double ky, double kz,
			     double x, double y, double z);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL*/
