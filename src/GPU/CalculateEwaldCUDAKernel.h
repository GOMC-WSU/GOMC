#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "XYZArray.h"

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
			  int imageSize);

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

__device__ double DotProduct(double kx, double ky, double kz,
			     double x, double y, double z);

#endif
