#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>

using namespace std;

__device__ int FlatIndexGPU(int i, int j);

void CallBoxInterForceGPU(vector<uint> pair1,
			  vector<uint> pair2,
			  XYZArray const &currentCoords,
			  XYZArray const &currentCOM,
			  BoxDimensions const& boxAxes,
			  MoleculeLookup const& molLookup,
			  Molecules const&mols,
			  bool electrostatic,
			  vector<double> particleCharge,
			  vector<int> particleKind,
			  vector<int> particleMol,
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
				 int *gpu_particleMol);

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
__device__ double CalcVirParticleGPU(double distSq, double index);
__device__ double CalcVirShiftGPU(double distSq, double index);
__device__ double CalcVirSwitchMartiniGPU(double distSq, double index);
__device__ double CalcVirSwitchGPU(double distSq, double index);

#endif
