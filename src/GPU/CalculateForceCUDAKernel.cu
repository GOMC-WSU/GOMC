#ifdef GOMC_CUDA

#include <cuda.h>
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "cub/cub.cuh"

using namespace cub;

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
			  uint const box)
{
  int atomNumber = currentCoords.Count();
  int molNumber = currentCOM.Count();
  int *gpu_pair1, *gpu_pair2;
  int *gpu_particleKind;
  int *gpu_particleMol;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_particleCharge;
  double *gpu_x, *gpu_y, *gpu_z;
  double *gpu_comx, *gpu_comy, *gpu_comz;
  double *gpu_virInter, *gpu_virReal;
  double *gpu_final_virInter, *gpu_final_virReal;

  CubDebugExit(cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int)));
  CubDebugExit(cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int)));
  CubDebugExit(cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_particleCharge, 
			  particleCharge.size() * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_particleKind, 
			  particleKind.size() * sizeof(int)));
  CubDebugExit(cudaMalloc((void**) &gpu_particleMol, 
			  particleMol.size() * sizeof(int)));
  CubDebugExit(cudaMalloc((void**) &gpu_virInter, 
			  pair1.size() * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_virReal, 
			  pair1.size() * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_final_virReal, 
			  pair1.size() * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_final_virInter, 
			  pair1.size() * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_comx, molNumber * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_comy, molNumber * sizeof(double)));
  CubDebugExit(cudaMalloc((void**) &gpu_comz, molNumber * sizeof(double)));

  CubDebugExit(cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_x, currentCoords.x, atomNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_y, currentCoords.y, atomNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_z, currentCoords.z, atomNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_comx, currentCOM.x, molNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_comy, currentCOM.y, molNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_comz, currentCOM.z, molNumber * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_particleCharge, &particleCharge[0],
			  particleCharge.size() * sizeof(double),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_particleKind, &particleKind[0],
			  particleKind.size() * sizeof(int),
			  cudaMemcpyHostToDevice));
  CubDebugExit(cudaMemcpy(gpu_particleMol, &particleMol[0],
			  particleMol.size() * sizeof(int),
			  cudaMemcpyHostToDevice));

  // Run the kernel...
  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size()/threadsPerBlock) + 1;
  BoxInterForceGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_pair1,
						       gpu_pair2,
						       gpu_x,
						       gpu_y,
						       gpu_z,
						       gpu_comx,
						       gpu_comy,
						       gpu_comz,
						       boxAxes.GetAxis(box).x,
						       boxAxes.GetAxis(box).y,
						       boxAxes.GetAxis(box).z,
						       electrostatic,
						       gpu_particleCharge,
						       gpu_particleKind,
						       gpu_particleMol,
						       gpu_virInter,
						       gpu_virReal,
						       pair1.size());

  // ReduceSum // Virial of LJ
  void * d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_virInter,
		    gpu_final_virInter, pair1.size());
  CubDebugExit(cudaMalloc(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_virInter,
		    gpu_final_virInter, pair1.size());
  cudaFree(d_temp_storage);

  // ReduceSum // Virial of Coulomb
  d_temp_storage = NULL;
  temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_virReal,
		    gpu_final_virReal, pair1.size());
  CubDebugExit(cudaMalloc(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_virReal,
		    gpu_final_virReal, pair1.size());
  cudaFree(d_temp_storage);
  
  // Copy back the result to CPU ! :)
  CubDebugExit(cudaMemcpy(&virInter, gpu_final_virInter, sizeof(double),
			  cudaMemcpyDeviceToHost));
  CubDebugExit(cudaMemcpy(&virReal, gpu_final_virReal, sizeof(double),
			  cudaMemcpyDeviceToHost));
  
  cudaFree(gpu_pair1);
  cudaFree(gpu_pair2);
  cudaFree(gpu_particleKind);
  cudaFree(gpu_particleMol);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  cudaFree(gpu_comx);
  cudaFree(gpu_comy);
  cudaFree(gpu_comz);
  cudaFree(gpu_virReal);
  cudaFree(gpu_virInter);
  cudaFree(gpu_final_virReal);
  cudaFree(gpu_final_virInter);
}

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
				 int pairSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID > pairSize)
    return;
  double distSq;
  double virX, virY, virZ;
  double pRF = 0.0, qi_qj, pVF = 0.0;
  //tensors for VDW and real part of electrostatic
  double vT11 = 0.0, vT12 = 0.0, vT13 = 0.0;
  double vT22 = 0.0, vT23 = 0.0, vT33 = 0.0;
  double rT11 = 0.0, rT12 = 0.0, rT13 = 0.0;
  double rT22 = 0.0, rT23 = 0.0, rT33 = 0.0;
  double diff_comx, diff_comy, diff_comz;

  if(InRcutGPU(distSq, virX, virY, virZ, gpu_x[gpu_pair1[threadID]], 
	       gpu_y[gpu_pair1[threadID]], gpu_z[gpu_pair1[threadID]], 
	       gpu_x[gpu_pair2[threadID]], gpu_y[gpu_pair2[threadID]], 
	       gpu_z[gpu_pair2[threadID]], xAxes, yAxes, zAxes, xAxes/2.0, 
	       yAxes/2.0, zAxes/2.0))
  {
    diff_comx = gpu_comx[gpu_particleMol[gpu_pair1[threadID]]] - 
      gpu_comx[gpu_particleMol[gpu_pair2[threadID]]];
    diff_comy = gpu_comy[gpu_particleMol[gpu_pair1[threadID]]] - 
      gpu_comy[gpu_particleMol[gpu_pair2[threadID]]];
    diff_comz = gpu_comz[gpu_particleMol[gpu_pair1[threadID]]] - 
      gpu_comz[gpu_particleMol[gpu_pair2[threadID]]];

    diff_comx = MinImageSignedGPU(diff_comx, xAxes, xAxes/2.0);
    diff_comy = MinImageSignedGPU(diff_comy, yAxes, yAxes/2.0);
    diff_comz = MinImageSignedGPU(diff_comz, zAxes, zAxes/2.0);

    if(electrostatic)
    {
      qi_qj = gpu_particleCharge[gpu_pair1[threadID]] * 
	gpu_particleCharge[gpu_pair2[threadID]];
      pRF = CalcCoulombForceGPU(distSq, qi_qj);
      
      rT11 = pRF * (virX * diff_comx);
      rT22 = pRF * (virY * diff_comy);
      rT33 = pRF * (virZ * diff_comz);
      
      //extra tensor calculations
      rT12 = pRF * (0.5 * (virX * diff_comy + virY * diff_comx));
      rT13 = pRF * (0.5 * (virX * diff_comz + virZ * diff_comx));
      rT23 = pRF * (0.5 * (virY * diff_comz + virZ * diff_comy));
    }
    
    pVF = CalcEnForceGPU(distSq, gpu_particleKind[gpu_pair1[threadID]], 
			 gpu_particleKind[gpu_pair2[threadID]]);
    
    vT11 = pVF * (virX * diff_comx);
    vT22 = pVF * (virY * diff_comy);
    vT33 = pVF * (virZ * diff_comz);
      
    //extra tensor calculations
    vT12 = pVF * (0.5 * (virX * diff_comy + virY * diff_comx));
    vT13 = pVF * (0.5 * (virX * diff_comz + virZ * diff_comx));
    vT23 = pVF * (0.5 * (virY * diff_comz + virZ * diff_comy));
  }

  gpu_virInter[threadID] = vT11 + vT22 + vT33;
  gpu_virReal[threadID] = rT11 + rT22 + rT33;
}

__device__ double CalcCoulombForceGPU(double distSq, double qi_qj)
{
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND)
  {
    return CalcCoulombVirParticleGPU(distSq, qi_qj);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND)
  {
    return CalcCoulombVirShiftGPU(distSq, qi_qj);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj);
  }
  else
    return CalcCoulombVirSwitchGPU(distSq, qi_qj);
}

__device__ double CalcEnForceGPU(double distSq, int kind1, int kind2)
{
  int index = FlatIndexGPU(kind1, kind2);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND)
  {
    return CalcVirParticleGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND)
  {
    return CalcVirShiftGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcVirSwitchMartiniGPU(distSq, index);
  }
  else
    return CalcVirSwitchGPU(distSq, index);
}

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj)
{
  double dist = sqrt(distSq);
  double constValue = 2.0 * gpu_alpha / sqrt(M_PI);
  double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
  double temp = 1.0 - erf(gpu_alpha * dist);
  return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
}

__device__ double CalcCoulombVirShiftGPU(double distSq, double qi_qj)
{
  if(gpu_ewald)
  {
    double dist = sqrt(distSq);
    double constValue = 2.0 * gpu_alpha / sqrt(M_PI);
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = 1.0 - erf(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  }
  else
  {
    double dist = sqrt(distSq);
    return qi_qj/(distSq * dist);
  }
}
__device__ double CalcCoulombVirSwitchMartiniGPU(double distSq, double qi_qj)
{
  if(gpu_ewald)
  {
     double dist = sqrt(distSq);
     double constValue = 2.0 * gpu_alpha / sqrt(M_PI);
     double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
     double temp = 1.0 - erf(gpu_alpha * dist);
     return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  }
  else
  {
     // in Martini, the Coulomb switching distance is zero, so we will have
     // sqrt(distSq) - rOnCoul =  sqrt(distSq)
     double dist = sqrt(distSq);
     double rij_ronCoul_2 = distSq;
     double rij_ronCoul_3 = dist * distSq;
     
     double A1 = 1.0 * (-(1.0+4)*gpu_rCut)/(pow(gpu_rCut,1.0+2) *
					   pow(gpu_rCut, 2));
     double B1 = -1.0 * (-(1.0+3)*gpu_rCut)/(pow(gpu_rCut,1.0+2) *
					    pow(gpu_rCut, 3));

     double virCoul = A1/rij_ronCoul_2 + B1/rij_ronCoul_3;
     return qi_qj * gpu_diElectric_1 * ( 1.0/(dist * distSq) + virCoul/dist);
  }
}

__device__ double CalcCoulombVirSwitchGPU(double distSq, double qi_qj)
{
  if(gpu_ewald)
  {
    double dist = sqrt(distSq);
    double constValue = 2.0 * gpu_alpha / sqrt(M_PI);
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = 1.0 - erf(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  }
  else
  {
    double rCutSq = gpu_rCut * gpu_rCut;
    double dist = sqrt(distSq);
    double switchVal = distSq/rCutSq - 1.0;
    switchVal *= switchVal;

    double dSwitchVal = 2.0 * (distSq/rCutSq - 1.0) * 2.0 * dist/rCutSq;
    return -1.0 * qi_qj * (dSwitchVal/distSq - switchVal/(distSq * dist));
  }
}

//VDW Calculation
//*****************************************************************//
__device__ double CalcVirParticleGPU(double distSq, int index)
{
  double rNeg2 = 1.0/distSq;
  double rRat2 = gpu_sigmaSq[index] * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  return gpu_epsilon_Cn[index] * 6.0 * 
    ((gpu_n[index]/6.0) * repulse - attract) * rNeg2;
}

__device__ double CalcVirShiftGPU(double distSq, int index)
{
  double rNeg2 = 1.0/distSq;
  double rRat2 = gpu_sigmaSq[index] * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  return gpu_epsilon_Cn[index] * 6.0 * 
    ((gpu_n[index]/6.0) * repulse - attract) * rNeg2;
}

__device__ double CalcVirSwitchMartiniGPU(double distSq, int index)
{
  double r_1 = 1.0/sqrt(distSq);
  double r_8 = pow(r_1, 8);
  double r_n2 = pow(r_1, gpu_n[index]+2);

  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;

  double pn = gpu_n[index];
  double An = pn * ((pn+1)*gpu_rOn - (pn+4)*gpu_rCut)/
    (pow(gpu_rCut, pn+2)*pow(gpu_rCut-gpu_rOn, 2));
  double Bn = -pn * ((pn+1)*gpu_rOn-(pn+3)*gpu_rCut)/
    (pow(gpu_rCut, pn+2)*pow(gpu_rCut-gpu_rOn, 3));

  double sig6 = pow(gpu_sigmaSq[index], 3);
  double sign = pow(gpu_sigmaSq[index], pn/2);

  double A6 = 6.0 * ((6.0+1)*gpu_rOn-(6.0+4)*gpu_rCut)/
    (pow(gpu_rCut,6.0+2)*pow(gpu_rCut-gpu_rOn, 2));
  double B6 = -6.0 * ((6.0+1)*gpu_rOn-(6.0+3)*gpu_rCut)/
    (pow(gpu_rCut,6.0+2)*pow(gpu_rCut-gpu_rOn, 3));

  double dshifttempRep = An * rij_ron_2 + Bn * rij_ron_3;
  double dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;
  
  const double dshiftRep = ( distSq > gpu_rOn * gpu_rOn ? 
			     dshifttempRep * r_1 : 0);
  const double dshiftAtt = ( distSq > gpu_rOn * gpu_rOn ?
			     dshifttempAtt * r_1 : 0);
  double Wij = gpu_epsilon_Cn[index] * (sign * (pn * r_n2 + dshiftRep) -
					sig6 * (6.0 * r_8 + dshiftAtt));
  return Wij;
}

__device__ double CalcVirSwitchGPU(double distSq, int index)
{
  double rCutSq = gpu_rCut * gpu_rCut;
  double rCutSq_rijSq = rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;
  double rOnSq = gpu_rOn * gpu_rOn;

  double rNeg2 = 1.0/distSq;
  double rRat2 = rNeg2 * gpu_sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  double factor1 = rCutSq - 3 * rOnSq;
  double factor2 = pow((rCutSq - rOnSq), -3);

  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
  double fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

  const double factE = ( distSq > rOnSq ? fE : 1.0);
  const double factW = ( distSq > rOnSq ? fW : 0.0);

  double Wij = gpu_epsilon_Cn[index] * 6.0 * 
    ((gpu_n[index]/6.0) * repulse - attract) * rNeg2;
  double Eij = gpu_epsilon_Cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

#endif
