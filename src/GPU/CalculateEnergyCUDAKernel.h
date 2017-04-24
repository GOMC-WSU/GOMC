#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"
#include "MoleculeLookup.h"

using namespace std;

void InitGPUForceField(double const *sigmaSq, double const *epsilon_Cn,
		       double const *n, uint VDW_Kind,
		       bool isMartini, int sizeSq);

__device__ double MinImageSignedGPU(double raw,double ax, double halfAx);

void CallBoxInterGPU(vector<uint> pair1,
		     vector<uint> pair2,
		     XYZArray const&coords,
		     BoxDimensions const& boxAxes,
		     MoleculeLookup const& molLookup,
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


__device__ bool InRcutGPU();

__device__ double CalcCoulombGPU();
__device__ double CalcEnGPU();


#endif
