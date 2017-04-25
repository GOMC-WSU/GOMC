#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"
#include "MoleculeLookup.h"

using namespace std;

__device__ int FlatIndexGPU(int i, int j);
__device__ double MinImageSignedGPU(double raw,double ax, double halfAx);

void CallBoxInterGPU(vector<uint> pair1,
		     vector<uint> pair2,
		     Coordinates const&coords,
		     BoxDimensions const& boxAxes,
		     MoleculeLookup const& molLookup,
		     Molecules const&mols,
		     bool electrostatic,
		     vector<double> particleCharge,
		     vector<int> particleKind,
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
			    int *gpu_particleKind);

// Call by calculate energy whether it is in rCut
__device__ bool InRcutGPU(double &distSq, double gpu_x1, double gpu_y1, 
			  double gpu_z1, double gpu_x2, double gpu_y2, 
			  double gpu_z2, double xAxes, double yAxes, 
			  double zAxes, double xHalfAxes, double yHalfAxes, 
			  double zHalfAxes);

// Call by force calculate to return the distance and virial component
__device__ bool InRcutGPU(double &distSq, double &virX, double &virY, 
			  double &virZ, double gpu_x1, double gpu_y1, 
			  double gpu_z1, double gpu_x2, double gpu_y2, 
			  double gpu_z2, double xAxes, double yAxes, 
			  double zAxes, double xHalfAxes, double yHalfAxes, 
			  double zHalfAxes);

__device__ double CalcCoulombGPU(double distSq, double qi_qj_fact);
__device__ double CalcEnGPU(double distSq, int kind1, int kind2);

__device__ double CalcEnParticleGPU(double distSq, int index);
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact);
__device__ double CalcEnShiftGPU(double distSq, int index);
__device__ double CalcCoulombShiftGPU(double distSq, double qi_qj_fact);
__device__ double CalcEnSwitchMartiniGPU(double distSq, int index);
__device__ double CalcCoulombSwitchMartiniGPU(double distSq, double qi_qj_fact);
__device__ double CalcEnSwitchGPU(double distSq, int index);
__device__ double CalcCoulombSwitchGPU(double distSq, double qi_qj_fact);

#endif
