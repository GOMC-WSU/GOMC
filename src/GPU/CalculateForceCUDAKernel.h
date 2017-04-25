#pragma once

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>

using namespace std;

void CallBoxInterForceGPU(vector<int> pair1,
			  vector<int> pair2,
			  Coordinates const& currentCoords,
			  COM const& currentCOM,
			  BoxDimensions const& boxAxes,
			  MoleculeLookup const& molLookup,
			  Molecules const&mols,
			  bool electrostatic,
			  vector<double> particleCharge,
			  vector<int> particleKind,
			  vector<int> particleMol,
			  uint const box);

#endif
