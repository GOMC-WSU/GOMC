#include "CalculateForceCUDAKernel.h"

#ifdef GOMC_CUDA

#include <cuda.h>
#include "ConstantDefinitionsCUDA.h"
#include "CalculateMinImageCUDA.h"

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
  BoxInterGPU<<<blocksPerGrid, threadsPerBlock>>>(gpu_pair1,
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
						  gpu_particleMol);

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
  delete [] cpu_x;
  delete [] cpu_y;
  delete [] cpu_z;
  delete [] cpu_comx;
  delete [] cpu_comy;
  delete [] cpu_comz;
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
				 int *gpu_particleMol)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  double distSq;
  double virX, virY, virZ;
  double pRF = 0.0, qi_qj, pVF = 0.0;
  //tensors for VDW and real part of electrostatic
  double vT11 = 0.0, vT12 = 0.0, vT13 = 0.0;
  double vT22 = 0.0, vT23 = 0.0, vT33 = 0.0;
  double rT11 = 0.0, rT12 = 0.0, rT13 = 0.0;
  double rT22 = 0.0, rT23 = 0.0, rT33 = 0.0;
  double diff_comx, diff_comy, diff_comz;

  if(InRcutGPU(&distSq, virX, virY, virZ, gpu_x[pair1[threadID]], 
	       gpu_y[pair1[threadID]], gpu_z[pair1[threadID]], 
	       gpu_x[pair2[threadID]], gpu_y[pair2[threadID]], 
	       gpu_z[pair2[threadID]], xAxes, yAxes, zAxes, xAxes/2.0, 
	       yAxes/2.0, zAxes/2.0))
  {
    diff_comx = gpu_comx[particleMol[pair1[threadID]]] - 
      gpu_comx[particleMol[pair2[threadID]]];
    diff_comy = gpu_comy[particleMol[pair1[threadID]]] - 
      gpu_comy[particleMol[pair2[threadID]]];
    diff_comz = gpu_comz[particleMol[pair1[threadID]]] - 
      gpu_comz[particleMol[pair2[threadID]]];

    diff_comx = MinImageSignedGPU(diff_comx, xAxes, xAxes/2.0);
    diff_comy = MinImageSignedGPU(diff_comy, yAxes, yAxes/2.0);
    diff_comz = MinImageSignedGPU(diff_comz, zAxes, zAxes/2.0);

    if(electrostatic)
    {
      qi_qj = gpu_particleCharge[gpu_pair1[threadID]] * 
	particleCharge[gpu_pair2[threadID]];
      pRF = CalcCoulombForceGPU(distSq, qi_qj);
      
      rT11 = pRF * (virX * diff_comx);
      rT22 = pRF * (virY * diff_comy);
      rT33 = pRF * (virZ * diff_comz);
      
      //extra tensor calculations
      rT12 = pRF * (0.5 * (virX * diff_comy + virY * diff_comx));
      rT13 = pRF * (0.5 * (virX * diff_comz + virZ * diff_comx));
      rT23 = pRF * (0.5 * (virY * diff_comz + virZ * diff_comy));
    }
    
    pVF = CalcEnForceGPU(distSq, particleKind[pair1[threadID]], 
			 particleKind[pair2[threadID]]);
    
    vT11 = pVF * (virX * diff_comx);
    vT22 = pVF * (virY * diff_comy);
    vT33 = pVF * (virZ * diff_comz);
      
    //extra tensor calculations
    vT12 = pVF * (0.5 * (virX * diff_comy + virY * diff_comx));
    vT13 = pVF * (0.5 * (virX * diff_comz + virZ * diff_comx));
    vT23 = pVF * (0.5 * (virY * diff_comz + virZ * diff_comy));
  }
}

__device__ double CalcCoulombForceGPU(double distSq, double qi_qj)
{
  if(gpu_VDW_Kind == VDW_STD_KIND)
  {
    return CalcCoulombVirParticleGPU(distSq, qi_qj);
  }
  else if(gpu_VDW_Kind == VDW_SHIFT_KIND)
  {
    return CalcCoulombVirShiftGPU(distSq, qi_qj);
  }
  else if(gpu_VDW_Kind == VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj);
  }
  else
    return CalcCoulombVirSwitchGPU(distSq, qi_qj);
}

__device__ double CalcEnForceGPU(double distSq, int kind1, int kind2)
{
  int index = FlatIndexGPU(kind1, kind2);
  if(gpu_VDW_Kind == VDW_STD_KIND)
  {
    return CalcVirParticleGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == VDW_SHIFT_KIND)
  {
    return CalcVirShiftGPU(distSq, index);
  }
  else if(gpu_VDW_Kind == VDW_SWITCH_KIND && gpu_isMartini)
  {
    return CalcVirSwitchMartiniGPU(distSq, index);
  }
  else
    return CalcVirSwitchGPU(distSq, index);
}

__device__ int FlatIndexGPU(int i, int j)
{
  return i + j * count;
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
  if(ewald)
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
     double rij_ronCoul_4 = distSq * distSq;
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
__device__ double CalcVirParticleGPU(double distSq, double index)
{
  double rNeg2 = 1.0/distSq;
  double rRat2 = gpu_sigmaSq[index] * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  return gpu_epsilon_Cn[index] * 6.0 * 
    ((gpu_n[index]/6.0) * repulse - attract) * rNeg2;
}

__device__ double CalcVirShiftGPU(double distSq, double index)
{
  double rNeg2 = 1.0/distSq;
  double rRat2 = gpu_sigmaSq[index] * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index]/2.0);
  return gpu_epsilon_Cn[index] * 6.0 * 
    ((gpu_n[index]/6.0) * repulse - attract) * rNeg2;
}

__device__ double CalcVirSwitchMartiniGPU(double distSq, double index)
{
  double r_1 = 1.0/sqrt(distSq);
  double r_8 = pow(r_1, 8);
  double r_n2 = pow(r1, gpu_n[index]+2);

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

__device__ double CalcVirSwitchGPU(double distSq, double index)
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
    (gpu_n[index]/6.0) * repulse - attract) * rNeg2;
  double Eij = gpu_epsilon_cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

#endif
