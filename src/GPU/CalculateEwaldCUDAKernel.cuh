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
			       double *sumInew,
			       double *prefact,
			       double &energyRecip);

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
			  double &energyRecipNew);

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
				      int imageSize);

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
				 int imageSize);

__global__ void BoxReciprocalGPU(double *gpu_prefact,
				 double *gpu_sumRnew,
				 double *gpu_sumInew,
				 double *gpu_energyRecip,
				 int imageSize);

__device__ double DotProduct(double kx, double ky, double kz,
			     double x, double y, double z);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_EWALD_CUDA_KERNEL*/
