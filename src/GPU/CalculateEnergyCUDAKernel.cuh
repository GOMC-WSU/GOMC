#ifndef CALCULATE_ENERGY_CUDA_KERNEL
#define CALCULATE_ENERGY_CUDA_KERNEL
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"

using namespace std;

void CallBoxInterGPU(vector<uint> pair1,
		     vector<uint> pair2,
		     XYZArray const &coords,
		     BoxDimensions const &boxAxes,
		     bool electrostatic,
		     vector<double> particleCharge,
		     vector<int> particleKind,
		     double &REn,
		     double &RLJEn,
		     uint const box);

__global__ void BoxInterGPU(int *gpu_pair1,
			    int *gpu_pair2,
			    double *gpu_x,
			    double *gpu_y,
			    double *gpu_z,
			    double xAxes,
			    double yAxes,
			    double zAxes,
			    bool electrostatic,
			    double *gpu_particleCharge,
			    int *gpu_particleKind,
			    double *gpu_REn,
			    double *gpu_LJEn,
			    int pairSize);


__device__ double CalcCoulombGPU(double distSq, double qi_qj_fact);
__device__ double CalcEnGPU(double distSq, int kind1, int kind2);

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact);
__device__ double CalcCoulombShiftGPU(double distSq, double qi_qj_fact);
__device__ double CalcCoulombSwitchMartiniGPU(double distSq, double qi_qj_fact);
__device__ double CalcCoulombSwitchGPU(double distSq, double qi_qj_fact);

//VDW Calculation
//*****************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index);
__device__ double CalcEnShiftGPU(double distSq, int index);
__device__ double CalcEnSwitchMartiniGPU(double distSq, int index);
__device__ double CalcEnSwitchGPU(double distSq, int index);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_ENERGY_CUDA_KERNEL*/
