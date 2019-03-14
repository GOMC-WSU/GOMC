/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA
#include "CalculateEnergyCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include <cuda.h>
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "cub/cub.cuh"
#include <stdio.h>

using namespace cub;

void CallBoxInterGPU(VariablesCUDA *vars,
                     vector<uint> pair1,
                     vector<uint> pair2,
                     XYZArray const &coords,
                     BoxDimensions const &boxAxes,
                     bool electrostatic,
                     vector<double> particleCharge,
                     vector<int> particleKind,
                     vector<int> particleMol,
                     double &REn,
                     double &LJEn,
                     bool multiParticleEnabled,
                     double *aForcex,
                     double *aForcey,
                     double *aForcez,
                     double *mForcex,
                     double *mForcey,
                     double *mForcez,
                     int atomCount,
                     int molCount,
                     bool reset_force,
                     bool copy_back,
                     uint const box)
{
  int atomNumber = coords.Count();
  int *gpu_pair1, *gpu_pair2, *gpu_particleKind, *gpu_particleMol;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_particleCharge;
  double *gpu_REn, *gpu_LJEn;
  double *gpu_final_REn, *gpu_final_LJEn;
  double cpu_final_REn, cpu_final_LJEn;

  if(reset_force) {
    cudaMemset(vars->gpu_aForcex, 0, atomCount * sizeof(double));
    cudaMemset(vars->gpu_aForcey, 0, atomCount * sizeof(double));
    cudaMemset(vars->gpu_aForcez, 0, atomCount * sizeof(double));
    cudaMemset(vars->gpu_mForcex, 0, molCount * sizeof(double));
    cudaMemset(vars->gpu_mForcey, 0, molCount * sizeof(double));
    cudaMemset(vars->gpu_mForcez, 0, molCount * sizeof(double));
  }

  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(double));
  cudaMalloc((void**) &gpu_particleKind, particleKind.size() * sizeof(int));
  cudaMalloc((void**) &gpu_particleMol, particleMol.size() * sizeof(int));
  cudaMalloc((void**) &gpu_REn, pair1.size() * sizeof(double));
  cudaMalloc((void**) &gpu_LJEn, pair1.size() * sizeof(double));
  cudaMalloc((void**) &gpu_final_REn, sizeof(double));
  cudaMalloc((void**) &gpu_final_LJEn, sizeof(double));
  
  // Copy necessary data to GPU
  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0],
             particleKind.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  // Run the kernel...
  threadsPerBlock = 256;
  blocksPerGrid = (int)(pair1.size() / threadsPerBlock) + 1;
  BoxInterGPU <<< blocksPerGrid, threadsPerBlock>>>(gpu_pair1,
      gpu_pair2,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      boxAxes.GetAxis(box).x,
      boxAxes.GetAxis(box).y,
      boxAxes.GetAxis(box).z,
      electrostatic,
      gpu_particleCharge,
      gpu_particleKind,
      gpu_particleMol,
      gpu_REn,
      gpu_LJEn,
      pair1.size(),
      vars->gpu_sigmaSq,
      vars->gpu_epsilon_Cn,
      vars->gpu_n,
      vars->gpu_VDW_Kind,
      vars->gpu_isMartini,
      vars->gpu_count,
      vars->gpu_rCut,
      vars->gpu_rCutCoulomb,
      vars->gpu_rCutLow,
      vars->gpu_rOn,
      vars->gpu_alpha,
      vars->gpu_ewald,
      vars->gpu_diElectric_1,
      vars->gpu_nonOrth,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      multiParticleEnabled,
      vars->gpu_aForcex,
      vars->gpu_aForcey,
      vars->gpu_aForcez,
      vars->gpu_mForcex,
      vars->gpu_mForcey,
      vars->gpu_mForcez,
      box);

  // ReduceSum
  void * d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                    gpu_final_REn, pair1.size());
  CubDebugExit(cudaMalloc(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                    gpu_final_REn, pair1.size());
  cudaFree(d_temp_storage);

  // LJ ReduceSum
  d_temp_storage = NULL;
  temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    gpu_final_LJEn, pair1.size());
  CubDebugExit(cudaMalloc(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    gpu_final_LJEn, pair1.size());
  cudaFree(d_temp_storage);
  // Copy back the result to CPU ! :)
  CubDebugExit(cudaMemcpy(&cpu_final_REn, gpu_final_REn, sizeof(double),
                          cudaMemcpyDeviceToHost));
  CubDebugExit(cudaMemcpy(&cpu_final_LJEn, gpu_final_LJEn, sizeof(double),
                          cudaMemcpyDeviceToHost));
  REn = cpu_final_REn;
  LJEn = cpu_final_LJEn;

  // Copy back forces
  if(multiParticleEnabled && copy_back) {
    CubDebugExit(cudaMemcpy(aForcex, vars->gpu_aForcex,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
    CubDebugExit(cudaMemcpy(aForcey, vars->gpu_aForcey,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
    CubDebugExit(cudaMemcpy(aForcez, vars->gpu_aForcez,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
    CubDebugExit(cudaMemcpy(mForcex, vars->gpu_mForcex,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
    CubDebugExit(cudaMemcpy(mForcey, vars->gpu_mForcey,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
    CubDebugExit(cudaMemcpy(mForcez, vars->gpu_mForcez,
                            sizeof(double) * atomCount,
                            cudaMemcpyDeviceToHost));
  }

  cudaFree(gpu_pair1);
  cudaFree(gpu_pair2);
  cudaFree(gpu_particleCharge);
  cudaFree(gpu_particleKind);
  cudaFree(gpu_REn);
  cudaFree(gpu_LJEn);
  cudaFree(gpu_final_REn);
  cudaFree(gpu_final_LJEn);
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
                            int *gpu_particleMol,
                            double *gpu_REn,
                            double *gpu_LJEn,
                            int pairSize,
                            double *gpu_sigmaSq,
                            double *gpu_epsilon_Cn,
                            double *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            double *gpu_rCut,
                            double *gpu_rCutCoulomb,
                            double *gpu_rCutLow,
                            double *gpu_rOn,
                            double *gpu_alpha,
                            int *gpu_ewald,
                            double *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            double *gpu_cell_x,
                            double *gpu_cell_y,
                            double *gpu_cell_z,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            bool multiParticleEnabled,
                            double *gpu_aForcex,
                            double *gpu_aForcey,
                            double *gpu_aForcez,
                            double *gpu_mForcex,
                            double *gpu_mForcey,
                            double *gpu_mForcez,
                            int box)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= pairSize)
    return;
  double distSq;
  double qi_qj_fact;
  double qqFact = 167000.0;
  double virX = 0.0, virY = 0.0, virZ = 0.0;
  double forceRealx = 0.0, forceRealy = 0.0, forceRealz = 0.0;
  double forceLJx = 0.0, forceLJy = 0.0, forceLJz = 0.0;
  gpu_REn[threadID] = 0.0;
  gpu_LJEn[threadID] = 0.0;
  double cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);
  if(InRcutGPU(distSq, virX, virY, virZ, gpu_x[gpu_pair1[threadID]],
               gpu_y[gpu_pair1[threadID]], gpu_z[gpu_pair1[threadID]],
               gpu_x[gpu_pair2[threadID]], gpu_y[gpu_pair2[threadID]],
               gpu_z[gpu_pair2[threadID]], xAxes, yAxes, zAxes, xAxes / 2.0,
               yAxes / 2.0, zAxes / 2.0, cutoff, gpu_nonOrth[0], gpu_cell_x,
               gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
               gpu_Invcell_z)) {
    if(electrostatic) {
      qi_qj_fact = gpu_particleCharge[gpu_pair1[threadID]] *
                   gpu_particleCharge[gpu_pair2[threadID]] * qqFact;
      gpu_REn[threadID] = CalcCoulombGPU(distSq, qi_qj_fact, gpu_rCutLow[0],
                                         gpu_ewald[0], gpu_VDW_Kind[0],
                                         gpu_alpha[box],
                                         gpu_rCutCoulomb[box],
                                         gpu_isMartini[0],
                                         gpu_diElectric_1[0]);
    }
    gpu_LJEn[threadID] = CalcEnGPU(distSq,
                                   gpu_particleKind[gpu_pair1[threadID]],
                                   gpu_particleKind[gpu_pair2[threadID]],
                                   gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                                   gpu_VDW_Kind[0], gpu_isMartini[0],
                                   gpu_rCut[0], gpu_rOn[0], gpu_count[0]);
    if(multiParticleEnabled) {
      if(electrostatic) {
        double coulombVir = CalcCoulombForceGPU(distSq, qi_qj_fact,
                                                gpu_VDW_Kind[0], gpu_ewald[0],
                                                gpu_isMartini[0],
                                                gpu_alpha[box],
                                                gpu_rCutCoulomb[box],
                                                gpu_diElectric_1[0]);
        forceRealx = virX * coulombVir;
        forceRealy = virY * coulombVir;
        forceRealz = virZ * coulombVir;
      }
      double pVF = CalcEnForceGPU(distSq, gpu_particleKind[gpu_pair1[threadID]],
                                  gpu_particleKind[gpu_pair2[threadID]],
                                  gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                                  gpu_rCut[0], gpu_rOn[0], gpu_isMartini[0],
                                  gpu_VDW_Kind[0], gpu_count[0]);
      forceLJx = virX * pVF;
      forceLJy = virY * pVF;
      forceLJz = virZ * pVF;

      atomicAdd(&gpu_aForcex[gpu_pair1[threadID]], forceRealx + forceLJx);
      atomicAdd(&gpu_aForcey[gpu_pair1[threadID]], forceRealy + forceLJy);
      atomicAdd(&gpu_aForcez[gpu_pair1[threadID]], forceRealz + forceLJz);
      atomicAdd(&gpu_aForcex[gpu_pair2[threadID]],-1.0*(forceRealx + forceLJx));
      atomicAdd(&gpu_aForcey[gpu_pair2[threadID]],-1.0*(forceRealy + forceLJy));
      atomicAdd(&gpu_aForcez[gpu_pair2[threadID]],-1.0*(forceRealz + forceLJz));

      atomicAdd(&gpu_mForcex[gpu_particleMol[gpu_pair1[threadID]]],
                forceRealx + forceLJx);
      atomicAdd(&gpu_mForcey[gpu_particleMol[gpu_pair1[threadID]]],
                forceRealy + forceLJy);
      atomicAdd(&gpu_mForcez[gpu_particleMol[gpu_pair1[threadID]]],
                forceRealz + forceLJz);
      atomicAdd(&gpu_mForcex[gpu_particleMol[gpu_pair2[threadID]]],
               -1.0 * (forceRealx + forceLJx));
      atomicAdd(&gpu_mForcey[gpu_particleMol[gpu_pair2[threadID]]],
               -1.0 * (forceRealy + forceLJy));
      atomicAdd(&gpu_mForcez[gpu_particleMol[gpu_pair2[threadID]]],
               -1.0 * (forceRealz + forceLJz));
    }
  }
}

__device__ double CalcCoulombGPU(double distSq, double qi_qj_fact,
                                 double gpu_rCutLow, int gpu_ewald,
                                 int gpu_VDW_Kind, double gpu_alpha,
                                 double gpu_rCutCoulomb, int gpu_isMartini,
                                 double gpu_diElectric_1)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }
  
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombParticleGPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcCoulombShiftGPU(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                               gpu_rCutCoulomb);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcCoulombSwitchMartiniGPU(distSq, qi_qj_fact, gpu_ewald,
                                       gpu_alpha, gpu_rCutCoulomb,
                                       gpu_diElectric_1);
  } else
    return CalcCoulombSwitchGPU(distSq, qi_qj_fact, gpu_alpha, gpu_ewald,
                                gpu_rCutCoulomb);
}

__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut, double gpu_rOn,
                            int gpu_count)
{
  if((gpu_rCut * gpu_rCut) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcEnParticleGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
  } else if(gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcEnShiftGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                          gpu_rCut);
  } else if(gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcEnSwitchMartiniGPU(distSq, index, gpu_sigmaSq, gpu_n,
                                  gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  } else
    return CalcEnSwitchGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                           gpu_rCut, gpu_rOn);
}

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, double qi_qj_fact,
                                         double gpu_ewald, double gpu_alpha)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * erfc(value) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_fact / dist;
  }
}

__device__ double CalcCoulombShiftGPU(double distSq, double qi_qj_fact,
                                      int gpu_ewald, double gpu_alpha,
                                      double gpu_rCut)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_fact * (1.0 / dist - 1.0 / gpu_rCut);
  }
}

__device__ double CalcCoulombSwitchMartiniGPU(double distSq, double qi_qj_fact,
                                              int gpu_ewald, double gpu_alpha,
                                              double gpu_rCut,
                                              double gpu_diElectric_1)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_3 = dist * distSq;
    double rij_ronCoul_4 = distSq * distSq;

    double A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 2));
    double B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 3));
    double C1 = 1.0 / pow(gpu_rCut, 1.0) - A1 / 3.0 * pow(gpu_rCut, 3) -
                B1 / 4.0 * pow(gpu_rCut, 4);

    double coul = -(A1 / 3.0) * rij_ronCoul_3 - (B1 / 4.0) * rij_ronCoul_4 - C1;
    return qi_qj_fact * gpu_diElectric_1 * (1.0 / dist + coul);
  }
}


__device__ double CalcCoulombSwitchGPU(double distSq, double qi_qj_fact,
                                       double gpu_alpha, int gpu_ewald,
                                       double gpu_rCut)
{
  if(gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    double rCutSq = gpu_rCut * gpu_rCut;
    double dist = sqrt(distSq);
    double switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_fact * switchVal / dist;
  }
}

//VDW Calculation
//**************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index,
                                    double *gpu_sigmaSq, double *gpu_n,
                                    double *gpu_epsilon_Cn)
{
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] / 2.0);
  return gpu_epsilon_Cn[index] * (repulse - attract);
}

__device__ double CalcEnShiftGPU(double distSq, int index, double *gpu_sigmaSq,
                                 double *gpu_n, double *gpu_epsilon_Cn,
                                 double gpu_rCut)
{
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] / 2.0);

  double shiftRRat2 = gpu_sigmaSq[index] / (gpu_rCut * gpu_rCut);
  double shiftRRat4 = shiftRRat2 * shiftRRat2;
  double shiftAttract = shiftRRat4 * shiftRRat2;
  double shiftRepulse = pow(shiftRRat2, gpu_n[index] / 2.0);
  double shiftConst = gpu_epsilon_Cn[index] * (shiftRepulse - shiftAttract);

  return (gpu_epsilon_Cn[index] * (repulse - attract) - shiftConst);
}

__device__ double CalcEnSwitchMartiniGPU(double distSq, int index,
    double *gpu_sigmaSq, double *gpu_n,
    double *gpu_epsilon_Cn,
    double gpu_rCut, double gpu_rOn)
{
  double r_2 = 1.0 / distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
  double r_n = pow(r_2, gpu_n[index] / 2.0);

  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double pn = gpu_n[index];
  double An = pn * ((pn + 1) * gpu_rOn - (pn + 4) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 2));
  double Bn = -pn * ((pn + 1) * gpu_rOn - (pn + 3) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 3));
  double Cn = 1.0 / pow(gpu_rCut, pn) - An / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              Bn / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  double A6 = 6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 4) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 2));
  double B6 = -6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 3) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 3));
  double C6 = 1.0 / pow(gpu_rCut, 6.0) - A6 / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              B6 / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  double shifttempRep = -(An / 3.0) * rij_ron_3 -
                        (Bn / 4.0) * rij_ron_4 - Cn;
  double shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const double shiftRep = ( distSq > gpu_rOn * gpu_rOn ? shifttempRep : -Cn);
  const double shiftAtt = ( distSq > gpu_rOn * gpu_rOn ? shifttempAtt : -C6);

  double sig6 = pow(gpu_sigmaSq[index], 3);
  double sign = pow(gpu_sigmaSq[index], pn / 2);
  double Eij = gpu_epsilon_Cn[index] * (sign * (r_n + shiftRep) -
                                        sig6 * (r_6 + shiftAtt));
  return Eij;
}


__device__ double CalcEnSwitchGPU(double distSq, int index, double *gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn,
                                  double gpu_rCut, double gpu_rOn)
{
  double rCutSq = gpu_rCut * gpu_rCut;
  double rOnSq = gpu_rOn * gpu_rOn;

  double rCutSq_rijSq = rCutSq  - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;

  double repulse = pow(rRat2, gpu_n[index] / 2.0);

  double factor1 = rCutSq - 3 * rOnSq;
  double factor2 = pow((rCutSq - rOnSq), -3);
  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const double factE = ( distSq > rOnSq ? fE : 1.0);

  return (gpu_epsilon_Cn[index] * (repulse - attract)) * factE;
}

#endif
