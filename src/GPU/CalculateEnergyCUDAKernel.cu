#include "CalculateEnergyCUDAKernel.h"

#ifdef GOMC_CUDA

#include <cuda.h>
__constant__ double gpu_sigmaSq[1000];
__constant__ double gpu_epsilon_Cn[1000];
__constant__ double gpu_n[1000];
__constant__ int gpu_VDW_Kind;
__constant__ bool gpu_isMartini;


void InitGPUForceField(double const *sigmaSq, double const *epsilon_Cn,
		       double const *n, uint VDW_Kind,
		       bool isMartini, int sizeSq)
{
  cudaMemcpyToSymbol("gpu_VDW_Kind", &VDW_Kind, sizeof(int));
  cudaMemcpyToSymbol("gpu_isMartini", &isMartini, sizeof(int));
  cudaMemcpyToSymbol("gpu_sigmaSq", sigmaSq, sizeSq * sizeof(double));
  cudaMemcpyToSymbol("gpu_epsilon_Cn", epsilon_Cn, sizeSq * sizeof(double));
  cudaMemcpyToSymbol("gpu_n", n, sizeSq * sizeof(double));
}

__device__ double MinImageSignedGpu(double raw,double ax, double halfAx) 
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

void CallBoxInterGPU(vector<int> pair1,
		     vector<int> pair2,
		     XYZArray const& coords,
		     BoxDimensions const& boxAxes,
		     MoleculeLookup const& molLookup,
		     bool electrostatic,
		     vector<double> particleCharge,
		     vector<int> particleKind,
		     uint const box)
{
  int molNumber = molLookup.NumInBox(box);
  int *gpu_pair1, *gpu_pair2;
  double *gpu_particleCharge;
  int *gpu_particleKind;
  double *gpu_x, *gpu_y, *gpu_z;
  double *cpu_x, *cpu_y, *cpu_z;
  int i=0;
  int blocksPerGrid, threadsPerBlock;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(box),
    end = molLookup.BoxEnd(box);

  cpu_x = new double[molNumber];
  cpu_y = new double[molNumber];
  cpu_z = new double[molNumber];
  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_x, molNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, molNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, molNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, particleCharge.size * sizeof(double));
  cudaMalloc((void**) &gpu_particleKind, particleKind.size * sizeof(int));

  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int), cudaMemcpyHosttoDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0], particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0], particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
  while(thisMol != end)
  {
    cpu_x[i] = coords.x[thisMol];
    cpu_y[i] = coords.y[thisMol];
    cpu_z[i] = coords.z[thisMol];
    i++;
    thisMol++;
  }
  cudaMemcpy(gpu_x, cpu_x, molNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, cpu_y, molNumber * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, cpu_z, molNumber * sizeof(double), cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size()/threadsPerBlock) + 1;
  BoxInterGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_pair1, gpu_pair2, gpu_x, gpu_y, gpu_z, boxAxes.GetAxis(box).x, boxAxes.GetAxis(box).y, boxAxes.GetAxis(box).z, electrostatic, gpu_particleCharge, gpu_particleKind);
  
  cudaFree(gpu_pair1);
  cudaFree(gpu_pair2);
  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  delete [] cpu_x;
  delete [] cpu_y;
  delete [] cpu_z;
}

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
			    int *gpu_particleKind)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  double distSq;
  double tempLJEn;
  double tmepREn;
  double qi_qj_fact;
  double qqFact = 167000.0;
  if(InRcutGPU(distSq, gpu_x, gpu_y, gpu_z, xAxes, yAxes, zAxes, xAxes/2.0, yAxes/2.0, zAxes/2.0))
  {
    if(electrostatic)
    {
      qi_qj_fact = gpu_particleCharge[gpu_pair1[threadID]] * particleCharge[gpu_pair2[threadID]] * qqFact;
      tempREn = CalcCoulombEnGPU(distSq, qi_qj_fact);
    }
    tempLJEn = CalcEnGPU(distSq, gpu_particleKind[gpu_pair1[threadID]], gpu_particleKind[gpu_pair2[threadID]]);
  } 
}

#endif
