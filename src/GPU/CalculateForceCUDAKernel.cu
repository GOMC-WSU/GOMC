#include "CalculateForceCUDAKernel.h"

#ifdef GOMC_CUDA

#include <cuda.h>
#include "ConstantDefinitionsCUDA.h"
#include "CalculateMinImageCUDA.h"

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
			  uint const box)
{
  int atomNumber = 0;
  for(int k = 0; k < molLookup.GetNumKind; k++)
  {
    atomNumber += molLookup.NumKindInBox(k, box) * mols.NumAtoms(k);
  }
  int molNumber = molLookup.NumInBox(box);
  int *gpu_pair1, *gpu_pair2;
  int start, length, i = 0, j = 0;
  int *gpu_particleKind;
  int *gpu_particleMol;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_particleCharge;
  double *gpu_x, *gpu_y, *gpu_z;
  double *cpu_x, *cpu_y, *cpu_z;
  double *gpu_comx, *gpu_comy, *gpu_comz;
  double *cpu_comx, *cpu_comy, *cpu_comz;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box),
    end = molLookup.BoxEnd(box);

  cpu_x = new double[atomNumber];
  cpu_y = new double[atomNumber];
  cpu_z = new double[atomNumber];
  cpu_comx = new double[molNumber];
  cpu_comy = new double[molNumber];
  cpu_comz = new double[molNumber];
  
  
  while(thisMol != end)
  {
    cpu_comx[i] = currentCOM.x[*thisMol];
    cpu_comy[i] = currentCOM.y[*thisMol];
    cpu_comz[i] = currentCOM.z[*thisMol];

    start = mols.MolStart(*thisMol);
    length = mols.NumAtomsByMol(*thisMol);
    
    for(int a = 0; a < length; a++)
    {
      cpu_x[j] = currentCoords.x[start + a];
      cpu_y[j] = currentCoords.y[start + a];
      cpu_z[j] = currentCoords.z[start + a];
      j++;
    }

    i++;
    thisMol++;
  }

  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size * sizeof(double));
  cudaMalloc((void**) &gpu_particleKind, particleKind.size * sizeof(int));
  cudaMalloc((void**) &gpu_particleMol, particleMol.size * sizeof(int));

  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int), 
	     cudaMemcpyHosttoDevice);
  cudaMemcpy(gpu_x, cpu_x, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, cpu_y, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, cpu_z, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_comx, cpu_comx, molNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_comy, cpu_comy, molNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_comz, cpu_comz, molNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0], 
	     particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0], 
	     particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], 
	     particleMol.size() * sizeof(int), cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size()/threadsPerBlock) + 1;

}

__global__ void BoxInterForceGPU()
{
}

#endif
