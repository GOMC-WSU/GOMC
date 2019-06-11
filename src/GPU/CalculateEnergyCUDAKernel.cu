/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA
#include "CalculateEnergyCUDAKernel.cuh"
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
                     vector<real> particleCharge,
                     vector<int> particleKind,
                     real &REn,
                     real &LJEn,
                     uint const box)
{
  int atomNumber = coords.Count();
  int *gpu_pair1, *gpu_pair2, *gpu_particleKind;
  int blocksPerGrid, threadsPerBlock;
  real *gpu_particleCharge;
  real *gpu_REn, *gpu_LJEn;
  real *gpu_final_REn, *gpu_final_LJEn;
  real cpu_final_REn, cpu_final_LJEn;

  cudaMalloc((void**) &gpu_pair1, pair1.size() * sizeof(int));
  cudaMalloc((void**) &gpu_pair2, pair2.size() * sizeof(int));
  cudaMalloc((void**) &gpu_particleCharge,
             particleCharge.size() * sizeof(real));
  cudaMalloc((void**) &gpu_particleKind,
             particleKind.size() * sizeof(int));
  cudaMalloc((void**) &gpu_REn, pair1.size() * sizeof(real));
  cudaMalloc((void**) &gpu_LJEn, pair1.size() * sizeof(real));
  cudaMalloc((void**) &gpu_final_REn, sizeof(real));
  cudaMalloc((void**) &gpu_final_LJEn, sizeof(real));


  // Copy necessary data to GPU
  cudaMemcpy(gpu_pair1, &pair1[0], pair1.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_pair2, &pair2[0], pair2.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0],
             particleKind.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(real),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(real),
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
  CubDebugExit(cudaMemcpy(&cpu_final_REn, gpu_final_REn, sizeof(real),
                          cudaMemcpyDeviceToHost));
  CubDebugExit(cudaMemcpy(&cpu_final_LJEn, gpu_final_LJEn, sizeof(real),
                          cudaMemcpyDeviceToHost));
  REn = cpu_final_REn;
  LJEn = cpu_final_LJEn;

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
                            real *gpu_x,
                            real *gpu_y,
                            real *gpu_z,
                            real xAxes,
                            real yAxes,
                            real zAxes,
                            bool electrostatic,
                            real *gpu_particleCharge,
                            int *gpu_particleKind,
                            real *gpu_REn,
                            real *gpu_LJEn,
                            int pairSize,
                            real *gpu_sigmaSq,
                            real *gpu_epsilon_Cn,
                            real *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            real *gpu_rCut,
                            real *gpu_rCutCoulomb,
                            real *gpu_rCutLow,
                            real *gpu_rOn,
                            real *gpu_alpha,
                            int *gpu_ewald,
                            real *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            real *gpu_cell_x,
                            real *gpu_cell_y,
                            real *gpu_cell_z,
                            real *gpu_Invcell_x,
                            real *gpu_Invcell_y,
                            real *gpu_Invcell_z,
                            int box)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= pairSize)
    return;
  real distSq;
  real qi_qj_fact;
  real qqFact = 167000.0;
  gpu_REn[threadID] = 0.0;
  gpu_LJEn[threadID] = 0.0;
  real cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);
  if(InRcutGPU(distSq, gpu_x[gpu_pair1[threadID]], gpu_y[gpu_pair1[threadID]],
               gpu_z[gpu_pair1[threadID]], gpu_x[gpu_pair2[threadID]],
               gpu_y[gpu_pair2[threadID]], gpu_z[gpu_pair2[threadID]],
               xAxes, yAxes, zAxes, xAxes / 2.0, yAxes / 2.0, zAxes / 2.0,
               cutoff, gpu_nonOrth[0], gpu_cell_x, gpu_cell_y,
               gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z)) {
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
  }
}

__device__ real CalcCoulombGPU(real distSq, real qi_qj_fact,
                                 real gpu_rCutLow, int gpu_ewald,
                                 int gpu_VDW_Kind, real gpu_alpha,
                                 real gpu_rCutCoulomb, int gpu_isMartini,
                                 real gpu_diElectric_1)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  if(gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombParticleGPU(distSq, qi_qj_fact, gpu_alpha);
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

__device__ real CalcEnGPU(real distSq, int kind1, int kind2,
                            real *gpu_sigmaSq, real *gpu_n,
                            real *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, real gpu_rCut, real gpu_rOn,
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
__device__ real CalcCoulombParticleGPU(real distSq, real qi_qj_fact,
    real gpu_alpha)
{
  real dist = sqrt(distSq);
  real value = gpu_alpha * dist;
  return qi_qj_fact * (1 - erf(value)) / dist;
}

__device__ real CalcCoulombShiftGPU(real distSq, real qi_qj_fact,
                                      int gpu_ewald, real gpu_alpha,
                                      real gpu_rCut)
{
  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    real dist = sqrt(distSq);
    return qi_qj_fact * (1.0 / dist - 1.0 / gpu_rCut);
  }
}

__device__ real CalcCoulombSwitchMartiniGPU(real distSq, real qi_qj_fact,
    int gpu_ewald, real gpu_alpha,
    real gpu_rCut,
    real gpu_diElectric_1)
{
  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    real dist = sqrt(distSq);
    real rij_ronCoul_3 = dist * distSq;
    real rij_ronCoul_4 = distSq * distSq;

    real A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 2));
    real B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
                pow(gpu_rCut, 3));
    real C1 = 1.0 / pow(gpu_rCut, 1.0) - A1 / 3.0 * pow(gpu_rCut, 3) -
                B1 / 4.0 * pow(gpu_rCut, 4);

    real coul = -(A1 / 3.0) * rij_ronCoul_3 - (B1 / 4.0) * rij_ronCoul_4 - C1;
    return qi_qj_fact * gpu_diElectric_1 * (1.0 / dist + coul);
  }
}


__device__ real CalcCoulombSwitchGPU(real distSq, real qi_qj_fact,
                                       real gpu_alpha, int gpu_ewald,
                                       real gpu_rCut)
{
  if(gpu_ewald) {
    real dist = sqrt(distSq);
    real value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    real rCutSq = gpu_rCut * gpu_rCut;
    real dist = sqrt(distSq);
    real switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_fact * switchVal / dist;
  }
}

//VDW Calculation
//**************************************************************//
__device__ real CalcEnParticleGPU(real distSq, int index,
                                    real *gpu_sigmaSq, real *gpu_n,
                                    real *gpu_epsilon_Cn)
{
  real rRat2 = gpu_sigmaSq[index] / distSq;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
  real repulse = pow(rRat2, gpu_n[index] / 2.0);
  return gpu_epsilon_Cn[index] * (repulse - attract);
}

__device__ real CalcEnShiftGPU(real distSq, int index, real *gpu_sigmaSq,
                                 real *gpu_n, real *gpu_epsilon_Cn,
                                 real gpu_rCut)
{
  real rRat2 = gpu_sigmaSq[index] / distSq;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
  real repulse = pow(rRat2, gpu_n[index] / 2.0);

  real shiftRRat2 = gpu_sigmaSq[index] / (gpu_rCut * gpu_rCut);
  real shiftRRat4 = shiftRRat2 * shiftRRat2;
  real shiftAttract = shiftRRat4 * shiftRRat2;
  real shiftRepulse = pow(shiftRRat2, gpu_n[index] / 2.0);
  real shiftConst = gpu_epsilon_Cn[index] * (shiftRepulse - shiftAttract);

  return (gpu_epsilon_Cn[index] * (repulse - attract) - shiftConst);
}

__device__ real CalcEnSwitchMartiniGPU(real distSq, int index,
    real *gpu_sigmaSq, real *gpu_n,
    real *gpu_epsilon_Cn,
    real gpu_rCut, real gpu_rOn)
{
  real r_2 = 1.0 / distSq;
  real r_4 = r_2 * r_2;
  real r_6 = r_4 * r_2;
  real r_n = pow(r_2, gpu_n[index] / 2.0);

  real rij_ron = sqrt(distSq) - gpu_rOn;
  real rij_ron_2 = rij_ron * rij_ron;
  real rij_ron_3 = rij_ron_2 * rij_ron;
  real rij_ron_4 = rij_ron_2 * rij_ron_2;

  real pn = gpu_n[index];
  real An = pn * ((pn + 1) * gpu_rOn - (pn + 4) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 2));
  real Bn = -pn * ((pn + 1) * gpu_rOn - (pn + 3) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 3));
  real Cn = 1.0 / pow(gpu_rCut, pn) - An / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              Bn / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  real A6 = 6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 4) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 2));
  real B6 = -6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 3) * gpu_rCut) /
              (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 3));
  real C6 = 1.0 / pow(gpu_rCut, 6.0) - A6 / 3.0 * pow(gpu_rCut - gpu_rOn, 3) -
              B6 / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  real shifttempRep = -(An / 3.0) * rij_ron_3 -
                        (Bn / 4.0) * rij_ron_4 - Cn;
  real shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const real shiftRep = ( distSq > gpu_rOn * gpu_rOn ? shifttempRep : -Cn);
  const real shiftAtt = ( distSq > gpu_rOn * gpu_rOn ? shifttempAtt : -C6);

  real sig6 = pow(gpu_sigmaSq[index], 3);
  real sign = pow(gpu_sigmaSq[index], pn / 2);
  real Eij = gpu_epsilon_Cn[index] * (sign * (r_n + shiftRep) -
                                        sig6 * (r_6 + shiftAtt));
  return Eij;
}


__device__ real CalcEnSwitchGPU(real distSq, int index, real *gpu_sigmaSq,
                                  real *gpu_n, real *gpu_epsilon_Cn,
                                  real gpu_rCut, real gpu_rOn)
{
  real rCutSq = gpu_rCut * gpu_rCut;
  real rOnSq = gpu_rOn * gpu_rOn;

  real rCutSq_rijSq = rCutSq  - distSq;
  real rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  real rRat2 = gpu_sigmaSq[index] / distSq;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;

  real repulse = pow(rRat2, gpu_n[index] / 2.0);

  real factor1 = rCutSq - 3 * rOnSq;
  real factor2 = pow((rCutSq - rOnSq), -3);
  real fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const real factE = ( distSq > rOnSq ? fE : 1.0);

  return (gpu_epsilon_Cn[index] * (repulse - attract)) * factE;
}

#endif
