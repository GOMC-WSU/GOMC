#ifndef CALCULATE_FORCE_CUDA_KERNEL
#define CALCULATE_FORCE_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"

using namespace std;

void CallBoxInterForceGPU(vector<uint> pair1,
			  vector<uint> pair2,
			  XYZArray const &currentCoords,
			  XYZArray const &currentCOM,
			  BoxDimensions const& boxAxes,
			  bool electrostatic,
			  vector<double> particleCharge,
			  vector<int> particleKind,
			  vector<int> particleMol,
			  double &virInter,
			  double &virReal,
			  uint const box);

__global__ void BoxInterForceGPU(int *gpu_pair1,
				 int *gpu_pair2,
				 double *gpu_x,
				 double *gpu_y,
				 double *gpu_z,
				 double *gpu_comx,
				 double *gpu_comy,
				 double *gpu_comz,
				 double xAxes,
				 double yAxes,
				 double zAxes,
				 bool electrostatic,
				 double *gpu_particleCharge,
				 int *gpu_particleKind,
				 int *gpu_particleMol,
				 double *gpu_virInter,
				 double *gpu_virReal,
				 int pairSize);

__device__ double CalcCoulombForceGPU(double distSq, double qi_qj);
__device__ double CalcEnForceGPU(double distSq, int kind1, int kind2);

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj);
__device__ double CalcCoulombVirShiftGPU(double distSq, double qi_qj);
__device__ double CalcCoulombVirSwitchMartiniGPU(double distSq, double qi_qj);
__device__ double CalcCoulombVirSwitchGPU(double distSq, double qi_qj);

//VDW Calculation
//*****************************************************************//
__device__ double CalcVirParticleGPU(double distSq, int index);
__device__ double CalcVirShiftGPU(double distSq, int index);
__device__ double CalcVirSwitchMartiniGPU(double distSq, int index);
__device__ double CalcVirSwitchGPU(double distSq, int index);

#endif /*GOMC_CUDA*/
#endif /*CALCULATE_FORCE_CUDA_KERNEL*/
