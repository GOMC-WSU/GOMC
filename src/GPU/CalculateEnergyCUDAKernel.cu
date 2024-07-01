/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA
#include <cuda.h>
#include <stdio.h>

#include <vector>

#include "CUDAMemoryManager.cuh"
#include "CalculateEnergyCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "cub/cub.cuh"
#define NUMBER_OF_NEIGHBOR_CELL 27

using namespace cub;

void CallBoxInterGPU(VariablesCUDA *vars, const std::vector<int> &cellVector,
                     const std::vector<int> &cellStartIndex,
                     const std::vector<std::vector<int>> &neighborList,
                     XYZArray const &coords, BoxDimensions const &boxAxes,
                     bool electrostatic,
                     const std::vector<double> &particleCharge,
                     const std::vector<int> &particleKind,
                     const std::vector<int> &particleMol, double &REn,
                     double &LJEn, bool sc_coul, double sc_sigma_6,
                     double sc_alpha, uint sc_power, uint const box) {
  int atomNumber = coords.Count();
  int neighborListCount = neighborList.size() * NUMBER_OF_NEIGHBOR_CELL;
  int numberOfCells = neighborList.size();
  int *gpu_particleKind, *gpu_particleMol;
  int *gpu_neighborList, *gpu_cellStartIndex;
  int blocksPerGrid, threadsPerBlock;
  int energyVectorLen;
  double *gpu_particleCharge;
  double *gpu_REn, *gpu_LJEn;

  // Run the kernel
  threadsPerBlock = 128;
  blocksPerGrid = numberOfCells * NUMBER_OF_NEIGHBOR_CELL;
  energyVectorLen = blocksPerGrid;

  // Convert neighbor list to 1D array
  std::vector<int> neighborlist1D(neighborListCount);
  for (int i = 0; i < neighborList.size(); i++) {
    for (int j = 0; j < NUMBER_OF_NEIGHBOR_CELL; j++) {
      neighborlist1D[i * NUMBER_OF_NEIGHBOR_CELL + j] = neighborList[i][j];
    }
  }

  CUMALLOC((void **)&gpu_neighborList, neighborListCount * sizeof(int));
  CUMALLOC((void **)&gpu_cellStartIndex, cellStartIndex.size() * sizeof(int));
  CUMALLOC((void **)&gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void **)&gpu_particleKind, particleKind.size() * sizeof(int));
  CUMALLOC((void **)&gpu_particleMol, particleMol.size() * sizeof(int));
  CUMALLOC((void **)&gpu_LJEn, energyVectorLen * sizeof(double));
  if (electrostatic) {
    CUMALLOC((void **)&gpu_REn, energyVectorLen * sizeof(double));
  }

  // Copy necessary data to GPU
  cudaMemcpy(gpu_neighborList, &neighborlist1D[0],
             neighborListCount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cellStartIndex, &cellStartIndex[0],
             cellStartIndex.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cellVector, &cellVector[0], atomNumber * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0],
             particleKind.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.GetAxis(box).x, boxAxes.GetAxis(box).y,
                              boxAxes.GetAxis(box).z);

  double3 halfAx =
      make_double3(boxAxes.GetAxis(box).x * 0.5, boxAxes.GetAxis(box).y * 0.5,
                   boxAxes.GetAxis(box).z * 0.5);

  BoxInterGPU<<<blocksPerGrid, threadsPerBlock>>>(
      gpu_cellStartIndex, vars->gpu_cellVector, gpu_neighborList, numberOfCells,
      vars->gpu_x, vars->gpu_y, vars->gpu_z, axis, halfAx, electrostatic,
      gpu_particleCharge, gpu_particleKind, gpu_particleMol, gpu_REn, gpu_LJEn,
      vars->gpu_sigmaSq, vars->gpu_epsilon_Cn, vars->gpu_n, vars->gpu_VDW_Kind,
      vars->gpu_isMartini, vars->gpu_count, vars->gpu_rCut,
      vars->gpu_rCutCoulomb, vars->gpu_rCutLow, vars->gpu_rOn, vars->gpu_alpha,
      vars->gpu_ewald, vars->gpu_diElectric_1, vars->gpu_nonOrth,
      vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box], sc_coul, sc_sigma_6, sc_alpha, sc_power,
      vars->gpu_rMin, vars->gpu_rMaxSq, vars->gpu_expConst, vars->gpu_molIndex,
      vars->gpu_lambdaVDW, vars->gpu_lambdaCoulomb, vars->gpu_isFraction, box);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  // LJ ReduceSum
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    vars->gpu_finalVal, energyVectorLen);
  CubDebugExit(CUMALLOC(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    vars->gpu_finalVal, energyVectorLen);
  // Copy back the result to CPU ! :)
  CubDebugExit(cudaMemcpy(&LJEn, vars->gpu_finalVal, sizeof(double),
                          cudaMemcpyDeviceToHost));
  if (electrostatic) {
    // Real Term ReduceSum
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                      vars->gpu_finalVal, energyVectorLen);
    // Copy back the result to CPU ! :)
    CubDebugExit(cudaMemcpy(&REn, vars->gpu_finalVal, sizeof(double),
                            cudaMemcpyDeviceToHost));
  } else {
    REn = 0.0;
  }
  
  CUFREE(d_temp_storage);
  CUFREE(gpu_particleCharge);
  CUFREE(gpu_particleKind);
  CUFREE(gpu_particleMol);
  CUFREE(gpu_LJEn);
  if (electrostatic) {
    CUFREE(gpu_REn);
  }
  CUFREE(gpu_neighborList);
  CUFREE(gpu_cellStartIndex);
}

__global__ void
BoxInterGPU(int *gpu_cellStartIndex, int *gpu_cellVector, int *gpu_neighborList,
            int numberOfCells, double *gpu_x, double *gpu_y, double *gpu_z,
            double3 axis, double3 halfAx, bool electrostatic,
            double *gpu_particleCharge, int *gpu_particleKind,
            int *gpu_particleMol, double *gpu_REn, double *gpu_LJEn,
            double *gpu_sigmaSq, double *gpu_epsilon_Cn, double *gpu_n,
            int *gpu_VDW_Kind, int *gpu_isMartini, int *gpu_count,
            double *gpu_rCut, double *gpu_rCutCoulomb, double *gpu_rCutLow,
            double *gpu_rOn, double *gpu_alpha, int *gpu_ewald,
            double *gpu_diElectric_1, int *gpu_nonOrth, double *gpu_cell_x,
            double *gpu_cell_y, double *gpu_cell_z, double *gpu_Invcell_x,
            double *gpu_Invcell_y, double *gpu_Invcell_z, bool sc_coul,
            double sc_sigma_6, double sc_alpha, uint sc_power, double *gpu_rMin,
            double *gpu_rMaxSq, double *gpu_expConst, int *gpu_molIndex,
            double *gpu_lambdaVDW, double *gpu_lambdaCoulomb,
            bool *gpu_isFraction, int box) {

  __shared__ double shr_cutoff;
  __shared__ int shr_particlesInsideCurrentCell, shr_numberOfPairs;
  __shared__ int shr_currentCellStartIndex, shr_neighborCellStartIndex;
  __shared__ bool shr_sameCell;

  int currentCell = blockIdx.x / NUMBER_OF_NEIGHBOR_CELL;
  int nCellIndex = blockIdx.x;
  int neighborCell = gpu_neighborList[nCellIndex];

  if (currentCell > neighborCell) {
    if (threadIdx.x == 0) {
	  gpu_LJEn[blockIdx.x] = 0.0;
      if (electrostatic) gpu_REn[blockIdx.x] = 0.0;
	}
    return;
  }
  
  if (threadIdx.x == 0) {
    // Calculate number of particles inside current Cell
    shr_currentCellStartIndex = gpu_cellStartIndex[currentCell];
    shr_particlesInsideCurrentCell =
        gpu_cellStartIndex[currentCell + 1] - shr_currentCellStartIndex;

    // Calculate number of particles inside neighbor Cell
    shr_neighborCellStartIndex = gpu_cellStartIndex[neighborCell];
    int particlesInsideNeighboringCell =
        gpu_cellStartIndex[neighborCell + 1] - shr_neighborCellStartIndex;

    shr_sameCell = currentCell == neighborCell;
    // Total number of pairs
    shr_numberOfPairs =
        shr_particlesInsideCurrentCell * particlesInsideNeighboringCell;
    shr_cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);
  }
  __syncthreads();

  // Specialize BlockReduce for a 1D block of 128 threads of type double
  using BlockReduce = cub::BlockReduce<double, 128>;

  // Allocate shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage LJEn_temp_storage;
  __shared__ typename BlockReduce::TempStorage REn_temp_storage;

  double LJEn = 0.0, REn = 0.0;

  for (int pairIndex = threadIdx.x; pairIndex < shr_numberOfPairs;
       pairIndex += blockDim.x) {
    int neighborParticleIndex = pairIndex / shr_particlesInsideCurrentCell;
    int currentParticleIndex = pairIndex % shr_particlesInsideCurrentCell;

    int currentParticle =
        gpu_cellVector[shr_currentCellStartIndex + currentParticleIndex];
    int neighborParticle = gpu_cellVector[shr_neighborCellStartIndex +
                                          neighborParticleIndex];

    int mA = gpu_particleMol[currentParticle];
    int mB = gpu_particleMol[neighborParticle];
    bool skip = mA == mB || (shr_sameCell && mA > mB);
    if (!skip) {
      // Check if they are within rcut
      double distSq = 0.0;

      if (InRcutGPU(distSq, gpu_x, gpu_y, gpu_z, currentParticle,
                    neighborParticle, axis, halfAx, shr_cutoff, gpu_nonOrth[0],
                    gpu_cell_x, gpu_cell_y, gpu_cell_z, gpu_Invcell_x,
                    gpu_Invcell_y, gpu_Invcell_z)) {
        int kA = gpu_particleKind[currentParticle];
        int kB = gpu_particleKind[neighborParticle];
 
        double lambdaVDW = DeviceGetLambdaVDW(mA, mB, box, gpu_isFraction,
                                              gpu_molIndex, gpu_lambdaVDW);
        LJEn += CalcEnGPU(
            distSq, kA, kB, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_VDW_Kind[0],
            gpu_isMartini[0], gpu_rCut[0], gpu_rOn[0], gpu_count[0], lambdaVDW,
            sc_sigma_6, sc_alpha, sc_power, gpu_rMin, gpu_rMaxSq, gpu_expConst);

        if (electrostatic) {
          double qi_qj_fact = gpu_particleCharge[currentParticle] *
                              gpu_particleCharge[neighborParticle];
          if (qi_qj_fact != 0.0) {
            qi_qj_fact *= qqFactGPU;
            double lambdaCoulomb = DeviceGetLambdaCoulomb(
                mA, mB, box, gpu_isFraction, gpu_molIndex, gpu_lambdaCoulomb);
            REn += CalcCoulombGPU(
                distSq, kA, kB, qi_qj_fact, gpu_rCutLow[0], gpu_ewald[0],
                gpu_VDW_Kind[0], gpu_alpha[box], gpu_rCutCoulomb[box],
                gpu_isMartini[0], gpu_diElectric_1[0], lambdaCoulomb, sc_coul,
                sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq, gpu_count[0]);
          }
        }
      }
    }
  }
  __syncthreads();

  // Compute the block-wide sum for thread 0
  double aggregate = BlockReduce(LJEn_temp_storage).Sum(LJEn);

  if (threadIdx.x == 0) {
    gpu_LJEn[blockIdx.x] = aggregate;
  }

  if (electrostatic) {
    // Need to sync the threads before reusing temp_storage
    // so using different variables
    // Compute the block-wide sum for thread 0
    aggregate = BlockReduce(REn_temp_storage).Sum(REn);
    if (threadIdx.x == 0)
      gpu_REn[blockIdx.x] = aggregate;
  }
}

__device__ double
CalcCoulombGPU(double distSq, int kind1, int kind2, double qi_qj_fact,
               double gpu_rCutLow, int gpu_ewald, int gpu_VDW_Kind,
               double gpu_alpha, double gpu_rCutCoulomb, int gpu_isMartini,
               double gpu_diElectric_1, double gpu_lambdaCoulomb, bool sc_coul,
               double sc_sigma_6, double sc_alpha, uint sc_power,
               double *gpu_sigmaSq, int gpu_count) {
  if ((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if (gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcCoulombParticleGPU(distSq, index, qi_qj_fact, gpu_ewald,
                                  gpu_alpha, gpu_lambdaCoulomb, sc_coul,
                                  sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq);
  } else if (gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcCoulombShiftGPU(distSq, index, qi_qj_fact, gpu_ewald, gpu_alpha,
                               gpu_rCutCoulomb, gpu_lambdaCoulomb, sc_coul,
                               sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq);
  } else if (gpu_VDW_Kind == GPU_VDW_EXP6_KIND) {
    return CalcCoulombExp6GPU(distSq, index, qi_qj_fact, gpu_ewald, gpu_alpha,
                              gpu_lambdaCoulomb, sc_coul, sc_sigma_6, sc_alpha,
                              sc_power, gpu_sigmaSq);
  } else if (gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcCoulombSwitchMartiniGPU(
        distSq, index, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCutCoulomb,
        gpu_diElectric_1, gpu_lambdaCoulomb, sc_coul, sc_sigma_6, sc_alpha,
        sc_power, gpu_sigmaSq);
  } else
    return CalcCoulombSwitchGPU(distSq, index, qi_qj_fact, gpu_alpha, gpu_ewald,
                                gpu_rCutCoulomb, gpu_lambdaCoulomb, sc_coul,
                                sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq);
}

__device__ double CalcEnGPU(double distSq, int kind1, int kind2,
                            double *gpu_sigmaSq, double *gpu_n,
                            double *gpu_epsilon_Cn, int gpu_VDW_Kind,
                            int gpu_isMartini, double gpu_rCut, double gpu_rOn,
                            int gpu_count, double gpu_lambdaVDW,
                            double sc_sigma_6, double sc_alpha, uint sc_power,
                            double *gpu_rMin, double *gpu_rMaxSq,
                            double *gpu_expConst) {
  if ((gpu_rCut * gpu_rCut) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if (gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcEnParticleGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                             gpu_lambdaVDW, sc_sigma_6, sc_alpha, sc_power);
  } else if (gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcEnShiftGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                          gpu_rCut, gpu_lambdaVDW, sc_sigma_6, sc_alpha,
                          sc_power);
  } else if (gpu_VDW_Kind == GPU_VDW_EXP6_KIND) {
    return CalcEnExp6GPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_lambdaVDW,
                         sc_sigma_6, sc_alpha, sc_power, gpu_rMin, gpu_rMaxSq,
                         gpu_expConst);
  } else if (gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcEnSwitchMartiniGPU(
        distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_rCut, gpu_rOn,
        gpu_lambdaVDW, sc_sigma_6, sc_alpha, sc_power);
  } else
    return CalcEnSwitchGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn,
                           gpu_rCut, gpu_rOn, gpu_lambdaVDW, sc_sigma_6,
                           sc_alpha, sc_power);
}

// ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombParticleGPU(double distSq, int index,
                                         double qi_qj_fact, int gpu_ewald,
                                         double gpu_alpha,
                                         double gpu_lambdaCoulomb, bool sc_coul,
                                         double sc_sigma_6, double sc_alpha,
                                         uint sc_power, double *gpu_sigmaSq) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombParticleGPUNoLambda(distSq, qi_qj_fact, gpu_ewald,
                                          gpu_alpha);
  }
  if (sc_coul) {
    double sigma6 =
        gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = cbrt(softDist6);
    return gpu_lambdaCoulomb * CalcCoulombParticleGPUNoLambda(
                                   softRsq, qi_qj_fact, gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombParticleGPUNoLambda(
                                   distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombParticleGPUNoLambda(double distSq,
                                                 double qi_qj_fact,
                                                 int gpu_ewald,
                                                 double gpu_alpha) {
  double dist = sqrt(distSq);
  double value = 1.0;
  if (gpu_ewald) {
    value = erfc(gpu_alpha * dist);
  }
  return qi_qj_fact * value / dist;
}

__device__ double CalcCoulombShiftGPU(double distSq, int index,
                                      double qi_qj_fact, int gpu_ewald,
                                      double gpu_alpha, double gpu_rCut,
                                      double gpu_lambdaCoulomb, bool sc_coul,
                                      double sc_sigma_6, double sc_alpha,
                                      uint sc_power, double *gpu_sigmaSq) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombShiftGPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                       gpu_rCut);
  }

  if (sc_coul) {
    double sigma6 =
        gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = cbrt(softDist6);
    return gpu_lambdaCoulomb * CalcCoulombShiftGPUNoLambda(softRsq, qi_qj_fact,
                                                           gpu_ewald, gpu_alpha,
                                                           gpu_rCut);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombShiftGPUNoLambda(distSq, qi_qj_fact,
                                                           gpu_ewald, gpu_alpha,
                                                           gpu_rCut);
  }
}

__device__ double CalcCoulombShiftGPUNoLambda(double distSq, double qi_qj_fact,
                                              int gpu_ewald, double gpu_alpha,
                                              double gpu_rCut) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1.0 - erf(value)) / dist;
  } else {
    return qi_qj_fact * (1.0 / dist - 1.0 / gpu_rCut);
  }
}

__device__ double CalcCoulombExp6GPU(double distSq, int index,
                                     double qi_qj_fact, int gpu_ewald,
                                     double gpu_alpha, double gpu_lambdaCoulomb,
                                     bool sc_coul, double sc_sigma_6,
                                     double sc_alpha, uint sc_power,
                                     double *gpu_sigmaSq) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombExp6GPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }

  if (sc_coul) {
    double sigma6 =
        gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = cbrt(softDist6);
    return gpu_lambdaCoulomb * CalcCoulombExp6GPUNoLambda(softRsq, qi_qj_fact,
                                                          gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombExp6GPUNoLambda(distSq, qi_qj_fact, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombExp6GPUNoLambda(double distSq, double qi_qj_fact,
                                             int gpu_ewald, double gpu_alpha) {
  double dist = sqrt(distSq);
  double value = 1.0;
  if (gpu_ewald) {
    value = erfc(gpu_alpha * dist);
  }
  return qi_qj_fact * value / dist;
}

__device__ double
CalcCoulombSwitchMartiniGPU(double distSq, int index, double qi_qj_fact,
                            int gpu_ewald, double gpu_alpha, double gpu_rCut,
                            double gpu_diElectric_1, double gpu_lambdaCoulomb,
                            bool sc_coul, double sc_sigma_6, double sc_alpha,
                            uint sc_power, double *gpu_sigmaSq) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombSwitchMartiniGPUNoLambda(
        distSq, qi_qj_fact, gpu_ewald, gpu_alpha, gpu_rCut, gpu_diElectric_1);
  }

  if (sc_coul) {
    double sigma6 =
        gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = cbrt(softDist6);
    return gpu_lambdaCoulomb * CalcCoulombSwitchMartiniGPUNoLambda(
                                   softRsq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                   gpu_rCut, gpu_diElectric_1);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombSwitchMartiniGPUNoLambda(
                                   distSq, qi_qj_fact, gpu_ewald, gpu_alpha,
                                   gpu_rCut, gpu_diElectric_1);
  }
}

__device__ double
CalcCoulombSwitchMartiniGPUNoLambda(double distSq, double qi_qj_fact,
                                    int gpu_ewald, double gpu_alpha,
                                    double gpu_rCut, double gpu_diElectric_1) {
  if (gpu_ewald) {
    double dist = sqrt(distSq);
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1 - erf(value)) / dist;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_3 = dist * distSq;
    double rij_ronCoul_4 = distSq * distSq;

    // Unoptimized computation
    // double A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
    // pow(gpu_rCut, 2));
    // double B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
    // pow(gpu_rCut, 3));
    // double C1 = 1.0 / pow(gpu_rCut, 1.0) - A1 / 3.0 * pow(gpu_rCut, 3) -
    // B1 / 4.0 * pow(gpu_rCut, 4);

    // Optimized computation
    double A1 = -5.0 / (gpu_rCut * gpu_rCut * gpu_rCut * gpu_rCut);
    double B1 = 4.0 / (gpu_rCut * gpu_rCut * gpu_rCut * gpu_rCut * gpu_rCut);
    double C1 = 1.0 / gpu_rCut - A1 / 3.0 * gpu_rCut * gpu_rCut * gpu_rCut -
                B1 / 4.0 * gpu_rCut * gpu_rCut * gpu_rCut * gpu_rCut;

    double coul = -(A1 / 3.0) * rij_ronCoul_3 - (B1 / 4.0) * rij_ronCoul_4 - C1;
    return qi_qj_fact * gpu_diElectric_1 * (1.0 / dist + coul);
  }
}

__device__ double CalcCoulombSwitchGPU(double distSq, int index,
                                       double qi_qj_fact, double gpu_alpha,
                                       int gpu_ewald, double gpu_rCut,
                                       double gpu_lambdaCoulomb, bool sc_coul,
                                       double sc_sigma_6, double sc_alpha,
                                       uint sc_power, double *gpu_sigmaSq) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombSwitchGPUNoLambda(distSq, qi_qj_fact, gpu_ewald,
                                        gpu_alpha, gpu_rCut);
  }

  if (sc_coul) {
    double sigma6 =
        gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 * dist6;
    double softRsq = cbrt(softDist6);
    return gpu_lambdaCoulomb *
           CalcCoulombSwitchGPUNoLambda(softRsq, qi_qj_fact, gpu_ewald,
                                        gpu_alpha, gpu_rCut);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombSwitchGPUNoLambda(distSq, qi_qj_fact, gpu_ewald,
                                        gpu_alpha, gpu_rCut);
  }
}

__device__ double CalcCoulombSwitchGPUNoLambda(double distSq, double qi_qj_fact,
                                               int gpu_ewald, double gpu_alpha,
                                               double gpu_rCut) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    double value = gpu_alpha * dist;
    return qi_qj_fact * (1.0 - erf(value)) / dist;
  } else {
    double rCutSq = gpu_rCut * gpu_rCut;
    double switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_fact * switchVal / dist;
  }
}

// VDW Calculation
//**************************************************************//
__device__ double CalcEnParticleGPU(double distSq, int index,
                                    double *gpu_sigmaSq, double *gpu_n,
                                    double *gpu_epsilon_Cn,
                                    double gpu_lambdaVDW, double sc_sigma_6,
                                    double sc_alpha, uint sc_power) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcEnParticleGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                     gpu_epsilon_Cn);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  return gpu_lambdaVDW * CalcEnParticleGPUNoLambda(softRsq, index, gpu_sigmaSq,
                                                   gpu_n, gpu_epsilon_Cn);
}

__device__ double CalcEnParticleGPUNoLambda(double distSq, int index,
                                            double *gpu_sigmaSq, double *gpu_n,
                                            double *gpu_epsilon_Cn) {
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] * 0.5);
  return gpu_epsilon_Cn[index] * (repulse - attract);
}

__device__ double CalcEnShiftGPU(double distSq, int index, double *gpu_sigmaSq,
                                 double *gpu_n, double *gpu_epsilon_Cn,
                                 double gpu_rCut, double gpu_lambdaVDW,
                                 double sc_sigma_6, double sc_alpha,
                                 uint sc_power) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcEnShiftGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                  gpu_epsilon_Cn, gpu_rCut);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  return gpu_lambdaVDW * CalcEnShiftGPUNoLambda(softRsq, index, gpu_sigmaSq,
                                                gpu_n, gpu_epsilon_Cn,
                                                gpu_rCut);
}

__device__ double CalcEnShiftGPUNoLambda(double distSq, int index,
                                         double *gpu_sigmaSq, double *gpu_n,
                                         double *gpu_epsilon_Cn,
                                         double gpu_rCut) {
  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] * 0.5);

  double shiftRRat2 = gpu_sigmaSq[index] / (gpu_rCut * gpu_rCut);
  double shiftRRat4 = shiftRRat2 * shiftRRat2;
  double shiftAttract = shiftRRat4 * shiftRRat2;
  double shiftRepulse = pow(shiftRRat2, gpu_n[index] * 0.5);
  double shiftConst = gpu_epsilon_Cn[index] * (shiftRepulse - shiftAttract);

  return (gpu_epsilon_Cn[index] * (repulse - attract) - shiftConst);
}

__device__ double CalcEnExp6GPU(double distSq, int index, double *gpu_sigmaSq,
                                double *gpu_n, double gpu_lambdaVDW,
                                double sc_sigma_6, double sc_alpha,
                                uint sc_power, double *gpu_rMin,
                                double *gpu_rMaxSq, double *gpu_expConst) {
  if (distSq < gpu_rMaxSq[index]) {
    return DBL_MAX;
  }
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcEnExp6GPUNoLambda(distSq, index, gpu_n, gpu_rMin, gpu_expConst);
  }
  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  return gpu_lambdaVDW *
         CalcEnExp6GPUNoLambda(softRsq, index, gpu_n, gpu_rMin, gpu_expConst);
}

__device__ double CalcEnExp6GPUNoLambda(double distSq, int index, double *gpu_n,
                                        double *gpu_rMin,
                                        double *gpu_expConst) {
  double dist = sqrt(distSq);
  double rRat = gpu_rMin[index] / dist;
  double rRat2 = rRat * rRat;
  double attract = rRat2 * rRat2 * rRat2;

  uint alph_ij = gpu_n[index];
  double repulse = (6.0 / alph_ij) * exp(alph_ij * (1.0 - 1.0 / rRat));
  return gpu_expConst[index] * (repulse - attract);
}

__device__ double
CalcEnSwitchMartiniGPU(double distSq, int index, double *gpu_sigmaSq,
                       double *gpu_n, double *gpu_epsilon_Cn, double gpu_rCut,
                       double gpu_rOn, double gpu_lambdaVDW, double sc_sigma_6,
                       double sc_alpha, uint sc_power) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcEnSwitchMartiniGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                          gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  }

  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  return gpu_lambdaVDW *
         CalcEnSwitchMartiniGPUNoLambda(softRsq, index, gpu_sigmaSq, gpu_n,
                                        gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
}

__device__ double
CalcEnSwitchMartiniGPUNoLambda(double distSq, int index, double *gpu_sigmaSq,
                               double *gpu_n, double *gpu_epsilon_Cn,
                               double gpu_rCut, double gpu_rOn) {
  double r_2 = 1.0 / distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
  double r_n = pow(r_2, gpu_n[index] * 0.5);

  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double pn = gpu_n[index];

  // Original Unoptimized computation
  // double An = pn * ((pn + 1) * gpu_rOn - (pn + 4) * gpu_rCut) /
  // (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 2));
  // double Bn = -pn * ((pn + 1) * gpu_rOn - (pn + 3) * gpu_rCut) /
  // (pow(gpu_rCut, pn + 2) * pow(gpu_rCut - gpu_rOn, 3));
  // double Cn = 1.0 / pow(gpu_rCut, pn) - An / 3.0 * pow(gpu_rCut - gpu_rOn, 3)
  // - Bn / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  // double A6 = 6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 4) * gpu_rCut) /
  // (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 2));
  // double B6 = -6.0 * ((6.0 + 1) * gpu_rOn - (6.0 + 3) * gpu_rCut) /
  // (pow(gpu_rCut, 6.0 + 2) * pow(gpu_rCut - gpu_rOn, 3));
  // double C6 = 1.0 / pow(gpu_rCut, 6.0) - A6 / 3.0 * pow(gpu_rCut - gpu_rOn,
  // 3) - B6 / 4.0 * pow(gpu_rCut - gpu_rOn, 4);

  // Optimized computation
  double An =
      pn * ((pn + 1.0) * gpu_rOn - (pn + 4.0) * gpu_rCut) /
      (pow(gpu_rCut, pn + 2.0) * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double Bn = -pn * ((pn + 1.0) * gpu_rOn - (pn + 3.0) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2.0) * (gpu_rCut - gpu_rOn) *
               (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double Cn = 1.0 / pow(gpu_rCut, pn) -
              An / 3.0 * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn) *
                  (gpu_rCut - gpu_rOn) -
              Bn / 4.0 * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn) *
                  (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn);

  double A6 =
      6.0 * (7.0 * gpu_rOn - 10.0 * gpu_rCut) /
      (pow(gpu_rCut, 8.0) * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double B6 = -6.0 * (7.0 * gpu_rOn - 9.0 * gpu_rCut) /
              (pow(gpu_rCut, 8.0) * (gpu_rCut - gpu_rOn) *
               (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double C6 = pow(gpu_rCut, -6.0) -
              A6 / 3.0 * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn) *
                  (gpu_rCut - gpu_rOn) -
              B6 / 4.0 * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn) *
                  (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn);

  double shifttempRep = -(An / 3.0) * rij_ron_3 - (Bn / 4.0) * rij_ron_4 - Cn;
  double shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const double shiftRep = (distSq > gpu_rOn * gpu_rOn ? shifttempRep : -Cn);
  const double shiftAtt = (distSq > gpu_rOn * gpu_rOn ? shifttempAtt : -C6);

  double sig6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  double sign = pow(gpu_sigmaSq[index], pn * 0.5);
  double Eij = gpu_epsilon_Cn[index] *
               (sign * (r_n + shiftRep) - sig6 * (r_6 + shiftAtt));
  return Eij;
}

__device__ double CalcEnSwitchGPU(double distSq, int index, double *gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn,
                                  double gpu_rCut, double gpu_rOn,
                                  double gpu_lambdaVDW, double sc_sigma_6,
                                  double sc_alpha, uint sc_power) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcEnSwitchGPUNoLambda(distSq, index, gpu_sigmaSq, gpu_n,
                                   gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  }
  double sigma6 = gpu_sigmaSq[index] * gpu_sigmaSq[index] * gpu_sigmaSq[index];
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  return gpu_lambdaVDW * CalcEnSwitchGPUNoLambda(softRsq, index, gpu_sigmaSq,
                                                 gpu_n, gpu_epsilon_Cn,
                                                 gpu_rCut, gpu_rOn);
}

__device__ double CalcEnSwitchGPUNoLambda(double distSq, int index,
                                          double *gpu_sigmaSq, double *gpu_n,
                                          double *gpu_epsilon_Cn,
                                          double gpu_rCut, double gpu_rOn) {
  double rCutSq = gpu_rCut * gpu_rCut;
  double rOnSq = gpu_rOn * gpu_rOn;

  double rCutSq_rijSq = rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  double rRat2 = gpu_sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;

  double repulse = pow(rRat2, gpu_n[index] * 0.5);

  double factor1 = rCutSq - 3.0 * rOnSq;
  double factor2 = pow((rCutSq - rOnSq), -3.0);
  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2.0 * distSq);

  const double factE = (distSq > rOnSq ? fE : 1.0);

  return (gpu_epsilon_Cn[index] * (repulse - attract)) * factE;
}

#endif
