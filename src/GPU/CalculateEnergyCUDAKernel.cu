#include "CalculateEnergyCUDAKernel.h"

#ifdef GOMC_CUDA

#include <cuda.h>
#include "ConstantDefinitionsCUDA.h"
#include "CalculateMinImageCUDA.h"

void CallBoxInterGPU(vector<int> pair1,
		     vector<int> pair2,
		     XYZArray const &coords,
		     BoxDimensions const &boxAxes,
		     bool electrostatic,
		     vector<double> particleCharge,
		     vector<int> particleKind,
		     uint const box)
{
  int atomNumber = coords.Count();
  int *gpu_pair1, *gpu_pair2, *gpu_particleKind;
  int start, length;
  int i = 0;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_particleCharge;
  double *gpu_x, *gpu_y, *gpu_z; 

  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_x, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_y, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_z, atomNumber * sizeof(double));
  cudaMalloc((void**) &gpu_particleCharge, 
	     particleCharge.size * sizeof(double));
  cudaMalloc((void**) &gpu_particleKind, particleKind.size * sizeof(int));

  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int), 
	     cudaMemcpyHosttoDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0], 
	     particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0], 
	     particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(gpu_x, coords.x, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_y, coords.y, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_z, coords.z, atomNumber * sizeof(double),
	     cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size()/threadsPerBlock) + 1;
  BoxInterGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_pair1, 
						  gpu_pair2, 
						  gpu_x, 
						  gpu_y, 
						  gpu_z, 
						  boxAxes.GetAxis(box).x, 
						  boxAxes.GetAxis(box).y, 
						  boxAxes.GetAxis(box).z, 
						  electrostatic, 
						  gpu_particleCharge, 
						  gpu_particleKind,
						  pair1.size());
  
  cudaFree(gpu_pair1);
  cudaFree(gpu_pair2);
  cudaFree(gpu_x);
  cudaFree(gpu_y);
  cudaFree(gpu_z);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_particleKind);
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
			    int *gpu_particleKind,
			    int pairSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID>=pairSize)
    return;
  double distSq;
  double tempLJEn;
  double tmepREn;
  double qi_qj_fact;
  double qqFact = 167000.0;
  if(InRcutGPU(distSq, gpu_x[pair1[threadID]], gpu_y[pair1[threadID]], 
	       gpu_z[pair1[threadID]], gpu_x[pair2[threadID]], 
	       gpu_y[pair2[threadID]], gpu_z[pair2[threadID]], 
	       xAxes, yAxes, zAxes, xAxes/2.0, yAxes/2.0, zAxes/2.0))
  {
    if(electrostatic)
    {
      qi_qj_fact = gpu_particleCharge[gpu_pair1[threadID]] * 
	particleCharge[gpu_pair2[threadID]] * qqFact;
      tempREn = CalcCoulombEnGPU(distSq, qi_qj_fact);
    }
    tempLJEn = CalcEnGPU(distSq, gpu_particleKind[gpu_pair1[threadID]], 
			 gpu_particleKind[gpu_pair2[threadID]]);
  } 
}

__device__ double CalcCoulombGPU(double distSq, double qi_qj_fact)
{
  double rCutLowSq = gpu_rCutLow * gpu_rCutLow;
  if(distSq <= rCutLowSq)
    return DBL_MAX;

  if(gpu_VDW_Kind == GPU_VDW_STD_KIND)
  {
    return CalcCoulombParticleGPU(distSq, qi_qj_fact);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND)
  {
    return CalcCoulombShiftGPU(distSq, qi_qj_fact);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcCoulombSwitchMartiniGPU(distSq, qi_qj_fact);
  }
  else
    return CalcCoulombSwitchGPU(distSq, qi_qj_fact);
}

__device__ double CalcEnGPU(double distSq, int kind1, int kind2)
{
  int index = FlatIndexGPU(kind1, kind2);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND)
  {
    return CalcEnParticleGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND)
  {
    return CalcEnShiftGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcEnSwitchMartiniGPU(distSq, index);
  }
  else
    return CalcEnSwitchGPU(distSq, index);
}

__device__ int FlatIndexGPU(int i, int j)
{
  return i + j * count;
}

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact)
{
  double dist = sqrt(distSq);
  double value = gpu_alpha * dist;
  return qi_qj_fact * (1 - erf(value))/dist;
}

__device__ double CalcCoulombShiftGPU(double distSq, double qi_qj_fact)
{
  if(gpu_ewald)
  {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value))/dist;
  }
  else
  {
    double dist = sqrt(distSq);
    return qi_qj_fact * (1.0/dist - 1.0/gpu_rCut);
  }
}

__device__ double CalcCoulombSwitchMartiniGPU(double distSq, double qi_qj_fact)
{
  if(gpu_ewald)
  {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  }
  else
  {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_3 = dist * distSq;
    double rij_ronCoul_4 = distSq * distSq;
    
    double A1 = 1.0 * (-(1.0+4)*gpu_rCut)/(pow(gpu_rCut,1.0+2) *
					   pow(gpu_rCut, 2));
    double B1 = -1.0 * (-(1.0+3)*gpu_rCut)/(pow(gpu_rCut,1.0+2) *
					    pow(gpu_rCut, 3));
    double C1 = 1.0/pow(gpu_rCut, 1.0) - A1/3.0*pow(gpu_rCut,3)-
      B1/4.0 *pow(gpu_rCut,4);

    double coul = -(A1/3.0) * rij_ronCoul_3 - (B1/4.0) * rij_ronCoul_4 - C1;
    return qi_qj_fact * gpu_diElectric_1 * (1.0/dist + coul);
  }
}


__device__ double CalcCoulombSwitchGPU(double distSq, double qi_qj_fact)
{
  if(gpu_ewald)
  {
    double dist = sqrt(distSq);
    double value = alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  }
  else
  {
    double rCutSq = gpu_rCut * gpu_rCut;
    double dist = sqrt(distSq);
    double switchVal = distSq/rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_fact * switchVal/dist;
  }
}

//VDW Calculation
//**************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index)
{
  double rRat2 = gpu_sigmaSq[index]/distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  return gpu_epsilon_Cn[index] * (repulse-attract);
}

__device__ double CalcEnShiftGPU(double distSq, int index)
{
  double rRat2 = gpu_sigmaSq[index]/distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);

  double shiftRRat2 = gpu_sigmaSq[index]/(gpu_rCut * gpu_rCut);
  double shiftRRat4 = shiftRRat2 * shiftRRat2;
  double shiftAttract = shiftRRat4 * shiftRRat2;
  double shiftRepulse = pow(shiftRRat2, gpu_n[index]/2.0);
  double shiftConst = gpu_epsilon_Cn[index] * (shiftRepulse - shiftAttract);

  return (gpu_epsilon_Cn[index] * (repulse-attract) - shiftConst);
}

__device__ double CalcEnSwitchMartiniGPU(double distSq, int index)
{
  double r_2 = 1.0/distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
  double r_n = pow(r_2, gpu_n[index]/2.0);
  
  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double pn = gpu_n[index];
  double An = pn * ((pn+1)*gpu_rOn - (pn+4)*gpu_rCut)/
    (pow(gpu_rCut, pn+2)*pow(gpu_rCut-gpu_rOn, 2));
  double Bn = -pn * ((pn+1)*gpu_rOn-(pn+3)*gpu_rCut)/
    (pow(gpu_rCut, pn+2)*pow(gpu_rCut-gpu_rOn, 3));
  double Cn = 1.0/pow(gpu_rCut, pn)-An/3.0*pow(gpu_rCut-gpu_rOn,3)-
    Bn/4.0*pow(gpu_rCut-gpu_rOn, 4);

  double A6 = 6.0 * ((6.0+1)*gpu_rOn-(6.0+4)*gpu_rCut)/
    (pow(gpu_rCut,6.0+2)*pow(gpu_rCut-gpu_rOn, 2));
  double B6 = -6.0 * ((6.0+1)*gpu_rOn-(6.0+3)*gpu_rCut)/
    (pow(gpu_rCut,6.0+2)*pow(gpu_rCut-gpu_rOn, 3));
  double C6 = 1.0/pow(gpu_rCut, 6.0)-A6/3.0*pow(gpu_rCut-gpu_rOn,3)-
    B6/4.0*pow(gpu_rCut-gpu_rOn,4);

  double shifttempRep = -(An/3.0)*rij_ron_3 -
    (Bn/4.0)*rij_ron_4 - Cn;
  double shifttempAtt = -(A6/3.0)*rij_ron_3 - (B6/4.0)*rij_ron_4 - C6;

  const double shiftRep = ( distSq > gpu_rOn*gpu_rOn ? shifttempRep : -Cn);
  const double shiftAtt = ( distSq > gpu_rOn*gpu_rOn ? shifttempAtt : -C6);

  double sig6 = pow(gpu_sigmaSq[index], 3);
  double sign = pow(gpu_sigmaSq[index], pn/2);
  double Eij = gpu_epsilon_Cn[index] * (sign * (r_n + shiftRep) - 
					sig6 * (r_6 + shiftAtt));
  return Eij;
}


__device__ double CalcEnSwitchGPU(double distSq, int index)
{
  double rCutSq = gpu_rCut * gpu_rCut;
  double rOnSq = gpu_rOn * gpu_rOn;
  
  double rCutSq_rijSq =rCutSq  - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rcutSq_rijSq;
  
  double rRat2 = gpu_sigmaSq[index]/distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  
  double repulse = power(rRat2, gpu_n[index]/2.0);
  
  double factor1 = rCutSq - 3 * rOnSq;
  double factor2 = pow((rCutSq - rOnSq), -3);
  double fE = rcutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const double facte = ( distSq > rOnSq ? fE : 1.0);
  
  return (gpu_epsilon_Cn[index] * (repulse-attract)) * factE;
}

#endif
