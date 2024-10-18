/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA

#include <cuda.h>

#include "CUDAMemoryManager.cuh"
#include "CalculateEnergyCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "cub/cub.cuh"

const int NUMBER_OF_NEIGHBOR_CELLS = 27;
const int PARTICLES_PER_BLOCK = 32;
const int THREADS_PER_BLOCK = 128;
const int THREADS_PER_BLOCK_SM = 64;

using namespace cub;

void CallBoxInterForceGPU(
    VariablesCUDA *vars, const std::vector<int> &cellVector,
    const std::vector<int> &cellStartIndex,
    const std::vector<std::vector<int>> &neighborList,
    const std::vector<int> &mapParticleToCell, XYZArray const &currentCoords,
    XYZArray const &currentCOM, BoxDimensions const &boxAxes,
    bool electrostatic, double &rT11, double &rT12, double &rT13, double &rT22,
    double &rT23, double &rT33, double &vT11, double &vT12, double &vT13,
    double &vT22, double &vT23, double &vT33, bool sc_coul, double sc_sigma_6,
    double sc_alpha, uint sc_power, uint const box) {
  int atomNumber = currentCoords.Count();
  int molNumber = currentCOM.Count();
  int numberOfCells = neighborList.size();
  int numberOfCellPairs = numberOfCells * NUMBER_OF_NEIGHBOR_CELLS;
  int *gpu_neighborList, *gpu_cellStartIndex;

  // Run the kernel...
  int threadsPerBlock = 128;
  int blocksPerGrid = numberOfCellPairs;
  int energyVectorLen = blocksPerGrid;

  // Convert neighbor list to 1D array
  std::vector<int> neighborlist1D(numberOfCellPairs);
  for (int i = 0; i < numberOfCells; i++) {
    for (int j = 0; j < NUMBER_OF_NEIGHBOR_CELLS; j++) {
      neighborlist1D[i * NUMBER_OF_NEIGHBOR_CELLS + j] = neighborList[i][j];
    }
  }

  CUMALLOC((void **)&gpu_neighborList, numberOfCellPairs * sizeof(int));
  CUMALLOC((void **)&gpu_cellStartIndex, cellStartIndex.size() * sizeof(int));
  UpdateEnergyVecs(vars, energyVectorLen, electrostatic);
#ifndef NDEBUG
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  cudaMemcpy(vars->gpu_mapParticleToCell, &mapParticleToCell[0],
             atomNumber * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_neighborList, &neighborlist1D[0],
             numberOfCellPairs * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cellStartIndex, &cellStartIndex[0],
             cellStartIndex.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cellVector, &cellVector[0], atomNumber * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, currentCOM.x, molNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, currentCOM.y, molNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, currentCOM.z, molNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.GetAxis(box).x, boxAxes.GetAxis(box).y,
                              boxAxes.GetAxis(box).z);

  double3 halfAx =
      make_double3(boxAxes.GetAxis(box).x * 0.5, boxAxes.GetAxis(box).y * 0.5,
                   boxAxes.GetAxis(box).z * 0.5);

  BoxInterForceGPU<<<blocksPerGrid, threadsPerBlock>>>(
      gpu_cellStartIndex, vars->gpu_cellVector, gpu_neighborList, numberOfCells,
      atomNumber, vars->gpu_mapParticleToCell, vars->gpu_x, vars->gpu_y,
      vars->gpu_z, vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, axis, halfAx,
      electrostatic, vars->gpu_particleCharge, vars->gpu_particleKind,
      vars->gpu_particleMol, vars->gpu_rT11, vars->gpu_rT12, vars->gpu_rT13,
      vars->gpu_rT22, vars->gpu_rT23, vars->gpu_rT33, vars->gpu_vT11,
      vars->gpu_vT12, vars->gpu_vT13, vars->gpu_vT22, vars->gpu_vT23,
      vars->gpu_vT33, vars->gpu_sigmaSq, vars->gpu_epsilon_Cn, vars->gpu_n,
      vars->gpu_VDW_Kind, vars->gpu_isMartini, vars->gpu_count, vars->gpu_rCut,
      vars->gpu_rCutCoulomb, vars->gpu_rCutLow, vars->gpu_rOn, vars->gpu_alpha,
      vars->gpu_ewald, vars->gpu_diElectric_1, vars->gpu_cell_x[box],
      vars->gpu_cell_y[box], vars->gpu_cell_z[box], vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box], vars->gpu_Invcell_z[box], vars->gpu_nonOrth,
      sc_coul, sc_sigma_6, sc_alpha, sc_power, vars->gpu_rMin, vars->gpu_rMaxSq,
      vars->gpu_expConst, vars->gpu_molIndex, vars->gpu_lambdaVDW,
      vars->gpu_lambdaCoulomb, vars->gpu_isFraction, box);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT11, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT11, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT12, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT12, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT13, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT13, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT22, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT22, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT23, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT23, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_vT33, vars->gpu_finalVal, energyVectorLen);
  cudaMemcpy(&vT33, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);

  if (electrostatic) {
    // ReduceSum // Virial of Coulomb
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT11, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT11, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT12, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT12, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT13, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT13, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT22, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT22, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT23, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT23, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
    DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                      vars->gpu_rT33, vars->gpu_finalVal, energyVectorLen);
    cudaMemcpy(&rT33, vars->gpu_finalVal, sizeof(double),
               cudaMemcpyDeviceToHost);
  }

  CUFREE(gpu_neighborList);
  CUFREE(gpu_cellStartIndex);
}

void CallBoxForceGPU(VariablesCUDA *vars, const std::vector<int> &cellVector,
                     const std::vector<int> &cellStartIndex,
                     const std::vector<std::vector<int>> &neighborList,
                     const std::vector<int> &mapParticleToCell,
                     XYZArray const &coords, BoxDimensions const &boxAxes,
                     bool electrostatic, double &REn, double &LJEn,
                     double *aForcex, double *aForcey, double *aForcez,
                     double *mForcex, double *mForcey, double *mForcez,
                     int atomCount, int molCount, bool sc_coul,
                     double sc_sigma_6, double sc_alpha, uint sc_power,
                     uint const box) {
  int atomNumber = coords.Count();
  int numberOfCells = neighborList.size();
  int numberOfCellPairs = numberOfCells * NUMBER_OF_NEIGHBOR_CELLS;
  int *gpu_neighborList, *gpu_cellStartIndex;
  double *gpu_REn, *gpu_LJEn;
  double cpu_final_REn = 0.0, cpu_final_LJEn = 0.0;

  int threadsPerBlock = THREADS_PER_BLOCK_SM;
  int blocksPerGrid = numberOfCells;
  int energyVectorLen = numberOfCells;

  // Convert neighbor list to 1D array
  std::vector<int> neighborlist1D(numberOfCellPairs);
  for (int i = 0; i < numberOfCells; i++) {
    for (int j = 0; j < NUMBER_OF_NEIGHBOR_CELLS; j++) {
      neighborlist1D[i * NUMBER_OF_NEIGHBOR_CELLS + j] = neighborList[i][j];
    }
  }

  CUMALLOC((void **)&gpu_neighborList, numberOfCellPairs * sizeof(int));
  CUMALLOC((void **)&gpu_cellStartIndex, cellStartIndex.size() * sizeof(int));
  CUMALLOC((void **)&gpu_LJEn, energyVectorLen * sizeof(double));
  if (electrostatic) {
    CUMALLOC((void **)&gpu_REn, energyVectorLen * sizeof(double));
  }

  // Copy necessary data to GPU
  cudaMemcpy(vars->gpu_aForcex, aForcex, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForcey, aForcey, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForcez, aForcez, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcex, mForcex, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcey, mForcey, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcez, mForcez, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mapParticleToCell, &mapParticleToCell[0],
             atomNumber * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_neighborList, &neighborlist1D[0],
             numberOfCellPairs * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_cellStartIndex, &cellStartIndex[0],
             cellStartIndex.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_cellVector, &cellVector[0], atomNumber * sizeof(int),
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

  BoxForceGPU<<<blocksPerGrid, threadsPerBlock>>>(
      gpu_cellStartIndex, vars->gpu_cellVector, gpu_neighborList, numberOfCells,
      atomNumber, vars->gpu_mapParticleToCell, vars->gpu_x, vars->gpu_y,
      vars->gpu_z, axis, halfAx, electrostatic, vars->gpu_particleCharge,
      vars->gpu_particleKind, vars->gpu_particleMol, gpu_REn, gpu_LJEn,
      vars->gpu_sigmaSq, vars->gpu_epsilon_Cn, vars->gpu_n, vars->gpu_VDW_Kind,
      vars->gpu_isMartini, vars->gpu_count, vars->gpu_rCut,
      vars->gpu_rCutCoulomb, vars->gpu_rCutLow, vars->gpu_rOn, vars->gpu_alpha,
      vars->gpu_ewald, vars->gpu_diElectric_1, vars->gpu_nonOrth,
      vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box], vars->gpu_aForcex, vars->gpu_aForcey,
      vars->gpu_aForcez, vars->gpu_mForcex, vars->gpu_mForcey,
      vars->gpu_mForcez, sc_coul, sc_sigma_6, sc_alpha, sc_power,
      vars->gpu_rMin, vars->gpu_rMaxSq, vars->gpu_expConst, vars->gpu_molIndex,
      vars->gpu_lambdaVDW, vars->gpu_lambdaCoulomb, vars->gpu_isFraction, box);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // LJ ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    vars->gpu_finalVal, energyVectorLen);
  CubDebugExit(CUMALLOC(&d_temp_storage, temp_storage_bytes));
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_LJEn,
                    vars->gpu_finalVal, energyVectorLen);
  // Copy the result back to CPU ! :)
  CubDebugExit(cudaMemcpy(&cpu_final_LJEn, vars->gpu_finalVal, sizeof(double),
                          cudaMemcpyDeviceToHost));
  if (electrostatic) {
    // Real Term ReduceSum
    DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_REn,
                      vars->gpu_finalVal, energyVectorLen);
    // Copy the result back to CPU ! :)
    CubDebugExit(cudaMemcpy(&cpu_final_REn, vars->gpu_finalVal, sizeof(double),
                            cudaMemcpyDeviceToHost));
  }
  CUFREE(d_temp_storage);

  REn = cpu_final_REn;
  LJEn = cpu_final_LJEn;

  cudaMemcpy(aForcex, vars->gpu_aForcex, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(aForcey, vars->gpu_aForcey, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(aForcez, vars->gpu_aForcez, sizeof(double) * atomCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(mForcex, vars->gpu_mForcex, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(mForcey, vars->gpu_mForcey, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(mForcez, vars->gpu_mForcez, sizeof(double) * molCount,
             cudaMemcpyDeviceToHost);
#ifndef NDEBUG
  cudaDeviceSynchronize();
#endif

  CUFREE(gpu_LJEn);
  if (electrostatic) {
    CUFREE(gpu_REn);
  }
  CUFREE(gpu_neighborList);
  CUFREE(gpu_cellStartIndex);
}

void CallVirialReciprocalGPU(VariablesCUDA *vars, XYZArray const &currentCoords,
                             XYZArray const &currentCOMDiff,
                             const std::vector<double> &molCharge, double &wT11,
                             double &wT12, double &wT13, double &wT22,
                             double &wT23, double &wT33, uint imageSize,
                             double constVal, uint box) {
  int atomNumber = currentCoords.Count();

  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dx, currentCOMDiff.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dy, currentCOMDiff.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_dz, currentCOMDiff.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_molCharge, &molCharge[0],
             molCharge.size() * sizeof(double), cudaMemcpyHostToDevice);

  // Initialize the virial terms to zero
  cudaMemset(vars->gpu_wT11, 0, imageSize * sizeof(double));
  cudaMemset(vars->gpu_wT12, 0, imageSize * sizeof(double));
  cudaMemset(vars->gpu_wT13, 0, imageSize * sizeof(double));
  cudaMemset(vars->gpu_wT22, 0, imageSize * sizeof(double));
  cudaMemset(vars->gpu_wT23, 0, imageSize * sizeof(double));
  cudaMemset(vars->gpu_wT33, 0, imageSize * sizeof(double));

  dim3 threadsPerBlock(128, 1, 1);
  int blocksPerGridX = (imageSize + threadsPerBlock.x - 1) / threadsPerBlock.x;
  int blocksPerGridY =
      (atomNumber + PARTICLES_PER_BLOCK - 1) / PARTICLES_PER_BLOCK;
  dim3 blocksPerGrid(blocksPerGridX, blocksPerGridY, 1);
  VirialReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_x, vars->gpu_y, vars->gpu_z, vars->gpu_dx, vars->gpu_dy,
      vars->gpu_dz, vars->gpu_kxRef[box], vars->gpu_kyRef[box],
      vars->gpu_kzRef[box], vars->gpu_prefactRef[box], vars->gpu_hsqrRef[box],
      vars->gpu_sumRref[box], vars->gpu_sumIref[box], vars->gpu_molCharge,
      vars->gpu_wT11, vars->gpu_wT12, vars->gpu_wT13, vars->gpu_wT22,
      vars->gpu_wT23, vars->gpu_wT33, constVal, imageSize, atomNumber);
#ifndef NDEBUG
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
#endif

  // ReduceSum -- Virial of Reciprocal
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT11, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT11, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT12, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT12, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT13, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT13, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT22, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT22, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT23, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT23, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
  DeviceReduce::Sum(vars->cub_reduce_storage, vars->cub_reduce_storage_size,
                    vars->gpu_wT33, vars->gpu_finalVal, imageSize);
  cudaMemcpy(&wT33, vars->gpu_finalVal, sizeof(double), cudaMemcpyDeviceToHost);
}

__global__ void BoxInterForceGPU(
    int *gpu_cellStartIndex, int *gpu_cellVector, int *gpu_neighborList,
    int numberOfCells, int atomNumber, int *gpu_mapParticleToCell,
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_comx,
    double *gpu_comy, double *gpu_comz, double3 axis, double3 halfAx,
    bool electrostatic, double *gpu_particleCharge, int *gpu_particleKind,
    int *gpu_particleMol, double *gpu_rT11, double *gpu_rT12, double *gpu_rT13,
    double *gpu_rT22, double *gpu_rT23, double *gpu_rT33, double *gpu_vT11,
    double *gpu_vT12, double *gpu_vT13, double *gpu_vT22, double *gpu_vT23,
    double *gpu_vT33, double *gpu_sigmaSq, double *gpu_epsilon_Cn,
    double *gpu_n, int *gpu_VDW_Kind, int *gpu_isMartini, int *gpu_count,
    double *gpu_rCut, double *gpu_rCutCoulomb, double *gpu_rCutLow,
    double *gpu_rOn, double *gpu_alpha, int *gpu_ewald,
    double *gpu_diElectric_1, double *gpu_cell_x, double *gpu_cell_y,
    double *gpu_cell_z, double *gpu_Invcell_x, double *gpu_Invcell_y,
    double *gpu_Invcell_z, int *gpu_nonOrth, bool sc_coul, double sc_sigma_6,
    double sc_alpha, uint sc_power, double *gpu_rMin, double *gpu_rMaxSq,
    double *gpu_expConst, int *gpu_molIndex, double *gpu_lambdaVDW,
    double *gpu_lambdaCoulomb, bool *gpu_isFraction, int box) {
  __shared__ double shr_cutoff;
  __shared__ int shr_particlesInsideCurrentCell, shr_numberOfPairs;
  __shared__ int shr_currentCellStartIndex, shr_neighborCellStartIndex;
  __shared__ bool shr_sameCell;
  double distSq;
  double3 virComponents;

  virComponents = make_double3(0.0, 0.0, 0.0);

  // tensors for VDW
  double local_vT11 = 0.0, local_vT22 = 0.0, local_vT33 = 0.0;
  double local_vT12 = 0.0, local_vT13 = 0.0, local_vT23 = 0.0;

  double local_rT11 = 0.0, local_rT22 = 0.0, local_rT33 = 0.0;
  double local_rT12 = 0.0, local_rT13 = 0.0, local_rT23 = 0.0;

  int currentCell = blockIdx.x / NUMBER_OF_NEIGHBOR_CELLS;
  int nCellIndex = blockIdx.x;
  int neighborCell = gpu_neighborList[nCellIndex];

  // Skip some cell pairs so we don't double count particle pairs
  if (currentCell > neighborCell) {
    if (threadIdx.x == 0) {
      gpu_vT11[blockIdx.x] = 0.0;
      gpu_vT22[blockIdx.x] = 0.0;
      gpu_vT33[blockIdx.x] = 0.0;
      gpu_vT12[blockIdx.x] = 0.0;
      gpu_vT13[blockIdx.x] = 0.0;
      gpu_vT23[blockIdx.x] = 0.0;

      if (electrostatic) {
        gpu_rT11[blockIdx.x] = 0.0;
        gpu_rT22[blockIdx.x] = 0.0;
        gpu_rT33[blockIdx.x] = 0.0;
        gpu_rT12[blockIdx.x] = 0.0;
        gpu_rT13[blockIdx.x] = 0.0;
        gpu_rT23[blockIdx.x] = 0.0;
      }
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

  double3 diff_com;

  for (int pairIndex = threadIdx.x; pairIndex < shr_numberOfPairs;
       pairIndex += blockDim.x) {
    int neighborParticleIndex = pairIndex / shr_particlesInsideCurrentCell;
    int currentParticleIndex = pairIndex % shr_particlesInsideCurrentCell;

    int currentParticle =
        gpu_cellVector[shr_currentCellStartIndex + currentParticleIndex];
    int neighborParticle =
        gpu_cellVector[shr_neighborCellStartIndex + neighborParticleIndex];

    // We don't process the same pair of cells twice, so we just need to check
    // to be sure we have different molecules. The exception is when the two
    // cells are the same, then we need to skip some pairs of molecules so we
    // don't double count any pairs.
    int mA = gpu_particleMol[currentParticle];
    int mB = gpu_particleMol[neighborParticle];
    bool skip = mA == mB || (shr_sameCell && mA > mB);
    if (!skip) {
      if (InRcutGPU(distSq, virComponents, gpu_x, gpu_y, gpu_z, currentParticle,
                    neighborParticle, axis, halfAx, shr_cutoff, gpu_nonOrth[0],
                    gpu_cell_x, gpu_cell_y, gpu_cell_z, gpu_Invcell_x,
                    gpu_Invcell_y, gpu_Invcell_z)) {
        int kA = gpu_particleKind[currentParticle];
        int kB = gpu_particleKind[neighborParticle];

        double lambdaVDW = DeviceGetLambdaVDW(mA, mB, box, gpu_isFraction,
                                              gpu_molIndex, gpu_lambdaVDW);

        diff_com = Difference3(gpu_comx, gpu_comy, gpu_comz, mA, mB);
        if (gpu_nonOrth[0])
          diff_com = MinImageNonOrthGPU(diff_com, axis, halfAx, gpu_cell_x,
                                        gpu_cell_y, gpu_cell_z, gpu_Invcell_x,
                                        gpu_Invcell_y, gpu_Invcell_z);
        else
          diff_com = MinImageGPU(diff_com, axis, halfAx);

        double pVF = CalcEnForceGPU(
            distSq, kA, kB, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_rCut[0],
            gpu_rOn[0], gpu_isMartini[0], gpu_VDW_Kind[0], gpu_count[0],
            lambdaVDW, sc_sigma_6, sc_alpha, sc_power, gpu_rMin, gpu_rMaxSq,
            gpu_expConst);

        local_vT11 += pVF * (virComponents.x * diff_com.x);
        local_vT22 += pVF * (virComponents.y * diff_com.y);
        local_vT33 += pVF * (virComponents.z * diff_com.z);

        // extra tensor calculations
        local_vT12 +=
            pVF * 0.5 *
            (virComponents.x * diff_com.y + virComponents.y * diff_com.x);
        local_vT13 +=
            pVF * 0.5 *
            (virComponents.x * diff_com.z + virComponents.z * diff_com.x);
        local_vT23 +=
            pVF * 0.5 *
            (virComponents.y * diff_com.z + virComponents.z * diff_com.y);

        if (electrostatic) {
          double qi_qj = gpu_particleCharge[currentParticle] *
                         gpu_particleCharge[neighborParticle];
          // skip particle pairs with no charge
          if (qi_qj != 0.0) {
            double lambdaCoulomb = DeviceGetLambdaCoulomb(
                mA, mB, box, gpu_isFraction, gpu_molIndex, gpu_lambdaCoulomb);
            double pRF = CalcCoulombForceGPU(
                distSq, qi_qj, gpu_VDW_Kind[0], gpu_ewald[0], gpu_isMartini[0],
                gpu_alpha[box], gpu_rCutCoulomb[box], gpu_diElectric_1[0],
                gpu_sigmaSq, sc_coul, sc_sigma_6, sc_alpha, sc_power,
                lambdaCoulomb, gpu_count[0], kA, kB);

            local_rT11 += pRF * (virComponents.x * diff_com.x);
            local_rT22 += pRF * (virComponents.y * diff_com.y);
            local_rT33 += pRF * (virComponents.z * diff_com.z);

            // extra tensor calculations
            local_rT12 +=
                pRF * 0.5 *
                (virComponents.x * diff_com.y + virComponents.y * diff_com.x);
            local_rT13 +=
                pRF * 0.5 *
                (virComponents.x * diff_com.z + virComponents.z * diff_com.x);
            local_rT23 +=
                pRF * 0.5 *
                (virComponents.y * diff_com.z + virComponents.z * diff_com.y);
          }
        }
      }
    }
  }
  __syncthreads();

  // Use BlockReduce to sum local tensor values across threads in the block
  using BlockReduce = cub::BlockReduce<double, THREADS_PER_BLOCK>;
  __shared__ typename BlockReduce::TempStorage vT11_temp_storage;
  __shared__ typename BlockReduce::TempStorage vT22_temp_storage;
  __shared__ typename BlockReduce::TempStorage vT33_temp_storage;
  __shared__ typename BlockReduce::TempStorage vT12_temp_storage;
  __shared__ typename BlockReduce::TempStorage vT13_temp_storage;
  __shared__ typename BlockReduce::TempStorage vT23_temp_storage;

  double aggregate_vT11 = BlockReduce(vT11_temp_storage).Sum(local_vT11);
  double aggregate_vT22 = BlockReduce(vT22_temp_storage).Sum(local_vT22);
  double aggregate_vT33 = BlockReduce(vT33_temp_storage).Sum(local_vT33);
  double aggregate_vT12 = BlockReduce(vT12_temp_storage).Sum(local_vT12);
  double aggregate_vT13 = BlockReduce(vT13_temp_storage).Sum(local_vT13);
  double aggregate_vT23 = BlockReduce(vT23_temp_storage).Sum(local_vT23);

  if (threadIdx.x == 0) {
    gpu_vT11[blockIdx.x] = aggregate_vT11;
    gpu_vT22[blockIdx.x] = aggregate_vT22;
    gpu_vT33[blockIdx.x] = aggregate_vT33;
    gpu_vT12[blockIdx.x] = aggregate_vT12;
    gpu_vT13[blockIdx.x] = aggregate_vT13;
    gpu_vT23[blockIdx.x] = aggregate_vT23;
  }
  if (electrostatic) {
    __shared__ typename BlockReduce::TempStorage rT11_temp_storage;
    __shared__ typename BlockReduce::TempStorage rT22_temp_storage;
    __shared__ typename BlockReduce::TempStorage rT33_temp_storage;
    __shared__ typename BlockReduce::TempStorage rT12_temp_storage;
    __shared__ typename BlockReduce::TempStorage rT13_temp_storage;
    __shared__ typename BlockReduce::TempStorage rT23_temp_storage;

    double aggregate_rT11 = BlockReduce(rT11_temp_storage).Sum(local_rT11);
    double aggregate_rT22 = BlockReduce(rT22_temp_storage).Sum(local_rT22);
    double aggregate_rT33 = BlockReduce(rT33_temp_storage).Sum(local_rT33);
    double aggregate_rT12 = BlockReduce(rT12_temp_storage).Sum(local_rT12);
    double aggregate_rT13 = BlockReduce(rT13_temp_storage).Sum(local_rT13);
    double aggregate_rT23 = BlockReduce(rT23_temp_storage).Sum(local_rT23);

    if (threadIdx.x == 0) {
      gpu_rT11[blockIdx.x] = aggregate_rT11;
      gpu_rT22[blockIdx.x] = aggregate_rT22;
      gpu_rT33[blockIdx.x] = aggregate_rT33;
      gpu_rT12[blockIdx.x] = aggregate_rT12;
      gpu_rT13[blockIdx.x] = aggregate_rT13;
      gpu_rT23[blockIdx.x] = aggregate_rT23;
    }
  }
}

__global__ void
BoxForceGPU(int *gpu_cellStartIndex, int *gpu_cellVector, int *gpu_neighborList,
            int numberOfCells, int atomNumber, int *gpu_mapParticleToCell,
            double *gpu_x, double *gpu_y, double *gpu_z, double3 axis,
            double3 halfAx, bool electrostatic, double *gpu_particleCharge,
            int *gpu_particleKind, int *gpu_particleMol, double *gpu_REn,
            double *gpu_LJEn, double *gpu_sigmaSq, double *gpu_epsilon_Cn,
            double *gpu_n, int *gpu_VDW_Kind, int *gpu_isMartini,
            int *gpu_count, double *gpu_rCut, double *gpu_rCutCoulomb,
            double *gpu_rCutLow, double *gpu_rOn, double *gpu_alpha,
            int *gpu_ewald, double *gpu_diElectric_1, int *gpu_nonOrth,
            double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
            double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
            double *gpu_aForcex, double *gpu_aForcey, double *gpu_aForcez,
            double *gpu_mForcex, double *gpu_mForcey, double *gpu_mForcez,
            bool sc_coul, double sc_sigma_6, double sc_alpha, uint sc_power,
            double *gpu_rMin, double *gpu_rMaxSq, double *gpu_expConst,
            int *gpu_molIndex, double *gpu_lambdaVDW, double *gpu_lambdaCoulomb,
            bool *gpu_isFraction, int box) {
  __shared__ double shr_cutoff;
  __shared__ int shr_particlesInsideCurrentCell, shr_currentCellStartIndex;
  double REn = 0.0, LJEn = 0.0;
  int currentCell = blockIdx.x;

  if (threadIdx.x == 0) {
    // Calculate number of particles inside current Cell
    shr_currentCellStartIndex = gpu_cellStartIndex[currentCell];
    shr_particlesInsideCurrentCell =
        gpu_cellStartIndex[currentCell + 1] - shr_currentCellStartIndex;
    shr_cutoff = fmax(gpu_rCut[0], gpu_rCutCoulomb[box]);
  }
  __syncthreads();

  for (int particleIdx = threadIdx.x; particleIdx < shr_particlesInsideCurrentCell; particleIdx += blockDim.x) {
    int particle = gpu_cellVector[shr_currentCellStartIndex + particleIdx];
    int mA = gpu_particleMol[particle];
    double3 forceComponents = make_double3(0.0, 0.0, 0.0);
    for (int neighborCellIdx = 0; neighborCellIdx < NUMBER_OF_NEIGHBOR_CELLS; ++neighborCellIdx) {
      int neighborCell = gpu_neighborList[currentCell * NUMBER_OF_NEIGHBOR_CELLS + neighborCellIdx];
      // Calculate number of particles inside neighbor cell
      int particlesInsideNeighboringCell =
          gpu_cellStartIndex[neighborCell + 1] - gpu_cellStartIndex[neighborCell];
      for (int neighborIdx = 0; neighborIdx < particlesInsideNeighboringCell; ++neighborIdx) {
        int neighbor = gpu_cellVector[gpu_cellStartIndex[neighborCell] + neighborIdx];
        int mB = gpu_particleMol[neighbor];
        // Check to be sure these are different molecules
        if (mA != mB) {
          double distSq;
          double3 virComponents;
          if (InRcutGPU(distSq, virComponents, gpu_x, gpu_y, gpu_z, particle,
                        neighbor, axis, halfAx, shr_cutoff, gpu_nonOrth[0],
                        gpu_cell_x, gpu_cell_y, gpu_cell_z, gpu_Invcell_x,
                        gpu_Invcell_y, gpu_Invcell_z)) {

            double lambdaVDW = DeviceGetLambdaVDW(mA, mB, box, gpu_isFraction,
                                                  gpu_molIndex, gpu_lambdaVDW);

            int kA = gpu_particleKind[particle];
            int kB = gpu_particleKind[neighbor];
            if (currentCell < neighborCell || (currentCell == neighborCell && particle < neighbor)) {
              LJEn += CalcEnGPU(
                  distSq, kA, kB, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_VDW_Kind[0],
                  gpu_isMartini[0], gpu_rCut[0], gpu_rOn[0], gpu_count[0], lambdaVDW,
                  sc_sigma_6, sc_alpha, sc_power, gpu_rMin, gpu_rMaxSq, gpu_expConst);
            }
            double forces = CalcEnForceGPU(
                distSq, kA, kB, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn, gpu_rCut[0],
                gpu_rOn[0], gpu_isMartini[0], gpu_VDW_Kind[0], gpu_count[0],
                lambdaVDW, sc_sigma_6, sc_alpha, sc_power, gpu_rMin, gpu_rMaxSq,
                gpu_expConst);
            double qi_qj_fact = 0.0;
            if (electrostatic) {
              qi_qj_fact = gpu_particleCharge[particle] *
                                  gpu_particleCharge[neighbor];
              if (qi_qj_fact != 0.0) {
                qi_qj_fact *= qqFactGPU;
                double lambdaCoulomb = DeviceGetLambdaCoulomb(
                    mA, mB, box, gpu_isFraction, gpu_molIndex, gpu_lambdaCoulomb);
                if (currentCell < neighborCell || (currentCell == neighborCell && particle < neighbor)) {
                    REn += CalcCoulombGPU(
                        distSq, kA, kB, qi_qj_fact, gpu_rCutLow[0], gpu_ewald[0],
                        gpu_VDW_Kind[0], gpu_alpha[box], gpu_rCutCoulomb[box],
                        gpu_isMartini[0], gpu_diElectric_1[0], lambdaCoulomb, sc_coul,
                        sc_sigma_6, sc_alpha, sc_power, gpu_sigmaSq, gpu_count[0]);
                }
                forces += CalcCoulombForceGPU(
                    distSq, qi_qj_fact, gpu_VDW_Kind[0], gpu_ewald[0],
                    gpu_isMartini[0], gpu_alpha[box], gpu_rCutCoulomb[box],
                    gpu_diElectric_1[0], gpu_sigmaSq, sc_coul, sc_sigma_6, sc_alpha,
                    sc_power, lambdaCoulomb, gpu_count[0], kA, kB);
              }
            }
            forceComponents.x += forces * virComponents.x;
            forceComponents.y += forces * virComponents.y;
            forceComponents.z += forces * virComponents.z;
          }
        }
      }
    }
    gpu_aForcex[particle] = forceComponents.x;
    gpu_aForcey[particle] = forceComponents.y;
    gpu_aForcez[particle] = forceComponents.z;

    atomicAdd(&gpu_mForcex[mA], forceComponents.x);
    atomicAdd(&gpu_mForcey[mA], forceComponents.y);
    atomicAdd(&gpu_mForcez[mA], forceComponents.z);
  }
  __syncthreads();

  // Specialize BlockReduce code
  using BlockReduce = cub::BlockReduce<double, THREADS_PER_BLOCK_SM>;

  // Allocating shared memory for BlockReduce
  __shared__ typename BlockReduce::TempStorage LJEn_temp_storage;

  // Compute the block-wide sum for thread 0
  double aggregate1 = BlockReduce(LJEn_temp_storage).Sum(LJEn);

  if (threadIdx.x == 0)
    gpu_LJEn[currentCell] = aggregate1;

  if (electrostatic) {
    // Need to sync the threads before reusing temp_storage
    // so using different variables
    __shared__ typename BlockReduce::TempStorage REn_temp_storage;

    // Compute the block-wide sum for thread 0
    double aggregate2 = BlockReduce(REn_temp_storage).Sum(REn);

    if (threadIdx.x == 0)
      gpu_REn[currentCell] = aggregate2;
  }
}

__global__ void VirialReciprocalGPU(
    double *gpu_x, double *gpu_y, double *gpu_z, double *gpu_comDx,
    double *gpu_comDy, double *gpu_comDz, double *gpu_kxRef, double *gpu_kyRef,
    double *gpu_kzRef, double *gpu_prefactRef, double *gpu_hsqrRef,
    double *gpu_sumRref, double *gpu_sumIref, double *gpu_molCharge,
    double *gpu_wT11, double *gpu_wT12, double *gpu_wT13, double *gpu_wT22,
    double *gpu_wT23, double *gpu_wT33, double constVal, uint imageSize,
    uint atomNumber) {
  __shared__ double shared_coords[PARTICLES_PER_BLOCK * 7];
  int imageID = blockIdx.x * blockDim.x + threadIdx.x;
  int offset_coordinates_index = blockIdx.y * PARTICLES_PER_BLOCK;
  int numberOfAtoms =
      min(PARTICLES_PER_BLOCK, atomNumber - offset_coordinates_index);
  if (threadIdx.x < numberOfAtoms) {
    shared_coords[threadIdx.x * 7] =
        gpu_x[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 1] =
        gpu_y[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 2] =
        gpu_z[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 3] =
        gpu_comDx[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 4] =
        gpu_comDy[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 5] =
        gpu_comDz[offset_coordinates_index + threadIdx.x];
    shared_coords[threadIdx.x * 7 + 6] =
        gpu_molCharge[offset_coordinates_index + threadIdx.x];
  }

  if (imageID >= imageSize)
    return;

  double wT11 = 0.0, wT12 = 0.0, wT13 = 0.0, wT22 = 0.0, wT23 = 0.0, wT33 = 0.0;
  double factor, dot;

  if (blockIdx.y == 0) {
    double constant_part = constVal + 1.0 / gpu_hsqrRef[imageID];
    factor =
        gpu_prefactRef[imageID] * (gpu_sumRref[imageID] * gpu_sumRref[imageID] +
                                   gpu_sumIref[imageID] * gpu_sumIref[imageID]);
    wT11 = factor * (1.0 - 2.0 * constant_part * gpu_kxRef[imageID] *
                               gpu_kxRef[imageID]);
    wT12 = factor *
           (-2.0 * constant_part * gpu_kxRef[imageID] * gpu_kyRef[imageID]);
    wT13 = factor *
           (-2.0 * constant_part * gpu_kxRef[imageID] * gpu_kzRef[imageID]);
    wT22 = factor * (1.0 - 2.0 * constant_part * gpu_kyRef[imageID] *
                               gpu_kyRef[imageID]);
    wT23 = factor *
           (-2.0 * constant_part * gpu_kyRef[imageID] * gpu_kzRef[imageID]);
    wT33 = factor * (1.0 - 2.0 * constant_part * gpu_kzRef[imageID] *
                               gpu_kzRef[imageID]);
  }
  __syncthreads();

  // Intramolecular part
#pragma unroll 32
  for (int particleID = 0; particleID < numberOfAtoms; particleID++) {
    double dotsin, dotcos;

    dot = DotProductGPU(gpu_kxRef[imageID], gpu_kyRef[imageID],
                        gpu_kzRef[imageID], shared_coords[particleID * 7],
                        shared_coords[particleID * 7 + 1],
                        shared_coords[particleID * 7 + 2]);
    sincos(dot, &dotsin, &dotcos);
    factor = gpu_prefactRef[imageID] * 2.0 * shared_coords[particleID * 7 + 6] *
             (gpu_sumIref[imageID] * dotcos - gpu_sumRref[imageID] * dotsin);

    wT11 += factor * (gpu_kxRef[imageID] * shared_coords[particleID * 7 + 3]);
    wT12 += factor * 0.5 *
            (gpu_kxRef[imageID] * shared_coords[particleID * 7 + 4] +
             gpu_kyRef[imageID] * shared_coords[particleID * 7 + 3]);
    wT13 += factor * 0.5 *
            (gpu_kxRef[imageID] * shared_coords[particleID * 7 + 5] +
             gpu_kzRef[imageID] * shared_coords[particleID * 7 + 3]);
    wT22 += factor * (gpu_kyRef[imageID] * shared_coords[particleID * 7 + 4]);
    wT23 += factor * 0.5 *
            (gpu_kyRef[imageID] * shared_coords[particleID * 7 + 5] +
             gpu_kzRef[imageID] * shared_coords[particleID * 7 + 4]);
    wT33 += factor * (gpu_kzRef[imageID] * shared_coords[particleID * 7 + 5]);
  }

  atomicAdd(&gpu_wT11[imageID], wT11);
  atomicAdd(&gpu_wT12[imageID], wT12);
  atomicAdd(&gpu_wT13[imageID], wT13);
  atomicAdd(&gpu_wT22[imageID], wT22);
  atomicAdd(&gpu_wT23[imageID], wT23);
  atomicAdd(&gpu_wT33[imageID], wT33);
}

__device__ double
CalcEnForceGPU(double distSq, int kind1, int kind2, double *gpu_sigmaSq,
               double *gpu_n, double *gpu_epsilon_Cn, double gpu_rCut,
               double gpu_rOn, int gpu_isMartini, int gpu_VDW_Kind,
               int gpu_count, double gpu_lambdaVDW, double sc_sigma_6,
               double sc_alpha, uint sc_power, double *gpu_rMin,
               double *gpu_rMaxSq, double *gpu_expConst) {
  if ((gpu_rCut * gpu_rCut) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  if (gpu_VDW_Kind == GPU_VDW_STD_KIND) {
    return CalcVirParticleGPU(distSq, index, gpu_sigmaSq[index], gpu_n,
                              gpu_epsilon_Cn, sc_sigma_6, sc_alpha, sc_power,
                              gpu_lambdaVDW);
  } else if (gpu_VDW_Kind == GPU_VDW_SHIFT_KIND) {
    return CalcVirShiftGPU(distSq, index, gpu_sigmaSq[index], gpu_n,
                           gpu_epsilon_Cn, sc_sigma_6, sc_alpha, sc_power,
                           gpu_lambdaVDW);
  } else if (gpu_VDW_Kind == GPU_VDW_EXP6_KIND) {
    return CalcVirExp6GPU(distSq, index, gpu_sigmaSq[index], gpu_n, gpu_rMin,
                          gpu_rMaxSq, gpu_expConst, sc_sigma_6, sc_alpha,
                          sc_power, gpu_lambdaVDW);
  } else if (gpu_VDW_Kind == GPU_VDW_SWITCH_KIND && gpu_isMartini) {
    return CalcVirSwitchMartiniGPU(
        distSq, index, gpu_sigmaSq[index], gpu_n, gpu_epsilon_Cn, gpu_rCut,
        gpu_rOn, sc_sigma_6, sc_alpha, sc_power, gpu_lambdaVDW);
  } else
    return CalcVirSwitchGPU(distSq, index, gpu_sigmaSq[index], gpu_epsilon_Cn,
                            gpu_n, gpu_rCut, gpu_rOn);
}

// ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj,
                                            int gpu_ewald, double gpu_alpha,
                                            int index, double gpu_sigmaSq,
                                            bool sc_coul, double sc_sigma_6,
                                            double sc_alpha, uint sc_power,
                                            double gpu_lambdaCoulomb) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombVirParticleGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }

  if (sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    return gpu_lambdaCoulomb * correction * correction *
           CalcCoulombVirParticleGPU(softRsq, qi_qj, gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombVirParticleGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombVirParticleGPU(const double distSq, const double qi_qj,
                                            const int gpu_ewald, const double gpu_alpha) {
  const double dist = sqrt(distSq);
  if (gpu_ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double alphaValue = gpu_alpha * M_2_SQRTPI;
    double expValue = exp(-gpu_alpha * gpu_alpha * distSq);
    double erfcValue = erfc(gpu_alpha * dist) / dist;
    return qi_qj * (alphaValue * expValue + erfcValue) / distSq;
  } else {
    return qi_qj / (distSq * dist);
  }
}

__device__ double CalcCoulombVirShiftGPU(double distSq, double qi_qj,
                                         int gpu_ewald, double gpu_alpha,
                                         int index, double gpu_sigmaSq,
                                         bool sc_coul, double sc_sigma_6,
                                         double sc_alpha, uint sc_power,
                                         double gpu_lambdaCoulomb) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombVirShiftGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }

  if (sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    return gpu_lambdaCoulomb * correction * correction *
           CalcCoulombVirShiftGPU(softRsq, qi_qj, gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombVirShiftGPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombVirShiftGPU(double distSq, double qi_qj,
                                         int gpu_ewald, double gpu_alpha) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double constValue = gpu_alpha * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = erfc(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    return qi_qj / (distSq * dist);
  }
}

__device__ double CalcCoulombVirExp6GPU(double distSq, double qi_qj,
                                        int gpu_ewald, double gpu_alpha,
                                        int index, double gpu_sigmaSq,
                                        bool sc_coul, double sc_sigma_6,
                                        double sc_alpha, uint sc_power,
                                        double gpu_lambdaCoulomb) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombVirExp6GPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }
  if (sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    return gpu_lambdaCoulomb * correction * correction *
           CalcCoulombVirExp6GPU(softRsq, qi_qj, gpu_ewald, gpu_alpha);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombVirExp6GPU(distSq, qi_qj, gpu_ewald, gpu_alpha);
  }
}

__device__ double CalcCoulombVirExp6GPU(double distSq, double qi_qj,
                                        int gpu_ewald, double gpu_alpha) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double constValue = gpu_alpha * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = erfc(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    return qi_qj / (distSq * dist);
  }
}

__device__ double CalcCoulombVirSwitchMartiniGPU(
    double distSq, double qi_qj, int gpu_ewald, double gpu_alpha,
    double gpu_rCut, double gpu_diElectric_1, int index, double gpu_sigmaSq,
    bool sc_coul, double sc_sigma_6, double sc_alpha, uint sc_power,
    double gpu_lambdaCoulomb) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                          gpu_rCut, gpu_diElectric_1);
  }

  if (sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    return gpu_lambdaCoulomb * correction * correction *
           CalcCoulombVirSwitchMartiniGPU(softRsq, qi_qj, gpu_ewald, gpu_alpha,
                                          gpu_rCut, gpu_diElectric_1);
  } else {
    return gpu_lambdaCoulomb *
           CalcCoulombVirSwitchMartiniGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                          gpu_rCut, gpu_diElectric_1);
  }
}

__device__ double CalcCoulombVirSwitchMartiniGPU(double distSq, double qi_qj,
                                                 int gpu_ewald,
                                                 double gpu_alpha,
                                                 double gpu_rCut,
                                                 double gpu_diElectric_1) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double constValue = gpu_alpha * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = erfc(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double rij_ronCoul_2 = 1.0 / distSq;
    double rij_ronCoul_3 = 1.0 / (dist * distSq);

    // Unoptimized version
    // double A1 = 1.0 * (-(1.0 + 4) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
    // pow(gpu_rCut, 2));
    // double B1 = -1.0 * (-(1.0 + 3) * gpu_rCut) / (pow(gpu_rCut, 1.0 + 2) *
    // pow(gpu_rCut, 3));
    // Optimized version
    double gpu_invrCut = 1.0 / gpu_rCut;
    double A1 = -5.0 * gpu_invrCut * gpu_invrCut * gpu_invrCut * gpu_invrCut;
    double B1 = 4.0 * gpu_invrCut * gpu_invrCut * gpu_invrCut * gpu_invrCut *
                gpu_invrCut;

    double virCoul = A1 * rij_ronCoul_2 + B1 * rij_ronCoul_3;
    return qi_qj * gpu_diElectric_1 * (rij_ronCoul_3 + virCoul / dist);
  }
}

__device__ double CalcCoulombVirSwitchGPU(double distSq, double qi_qj,
                                          int gpu_ewald, double gpu_alpha,
                                          double gpu_rCut, int index,
                                          double gpu_sigmaSq, bool sc_coul,
                                          double sc_sigma_6, double sc_alpha,
                                          uint sc_power,
                                          double gpu_lambdaCoulomb) {
  if (gpu_lambdaCoulomb >= 0.999999) {
    return CalcCoulombVirSwitchGPU(distSq, qi_qj, gpu_ewald, gpu_alpha,
                                   gpu_rCut);
  }

  if (sc_coul) {
    double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
    sigma6 = max(sigma6, sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        sc_alpha * pow((1.0 - gpu_lambdaCoulomb), (double)sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    return gpu_lambdaCoulomb * correction * correction *
           CalcCoulombVirSwitchGPU(softRsq, qi_qj, gpu_ewald, gpu_alpha,
                                   gpu_rCut);
  } else {
    return gpu_lambdaCoulomb * CalcCoulombVirSwitchGPU(distSq, qi_qj, gpu_ewald,
                                                       gpu_alpha, gpu_rCut);
  }
}

__device__ double CalcCoulombVirSwitchGPU(double distSq, double qi_qj,
                                          int gpu_ewald, double gpu_alpha,
                                          double gpu_rCut) {
  double dist = sqrt(distSq);
  if (gpu_ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double constValue = gpu_alpha * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * gpu_alpha * gpu_alpha * distSq);
    double temp = erfc(gpu_alpha * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double rCutSq = gpu_rCut * gpu_rCut;
    double switchVal = distSq / rCutSq - 1.0;
    switchVal *= switchVal;

    double dSwitchVal = 2.0 * (distSq / rCutSq - 1.0) * 2.0 * dist / rCutSq;
    return -1.0 * qi_qj * (dSwitchVal / distSq - switchVal / (distSq * dist));
  }
}

// VDW Calculation
//*****************************************************************//
__device__ double CalcVirParticleGPU(double distSq, int index,
                                     double gpu_sigmaSq, double *gpu_n,
                                     double *gpu_epsilon_Cn, double sc_sigma_6,
                                     double sc_alpha, uint sc_power,
                                     double gpu_lambdaVDW) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcVirParticleGPU(distSq, index, gpu_sigmaSq, gpu_n,
                              gpu_epsilon_Cn);
  }

  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  return gpu_lambdaVDW * correction * correction *
         CalcVirParticleGPU(softRsq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
}

__device__ double CalcVirParticleGPU(double distSq, int index,
                                     double gpu_sigmaSq, double *gpu_n,
                                     double *gpu_epsilon_Cn) {
  double rNeg2 = 1.0 / distSq;
  double rRat2 = gpu_sigmaSq * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] * 0.5);
  // return gpu_epsilon_Cn[index] * 6.0 *
  // ((gpu_n[index] / 6.0) * repulse - attract) * rNeg2;
  return gpu_epsilon_Cn[index] * rNeg2 *
         (gpu_n[index] * repulse - 6.0 * attract);
}

__device__ double CalcVirShiftGPU(double distSq, int index, double gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn,
                                  double sc_sigma_6, double sc_alpha,
                                  uint sc_power, double gpu_lambdaVDW) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcVirShiftGPU(distSq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
  }

  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  return gpu_lambdaVDW * correction * correction *
         CalcVirShiftGPU(softRsq, index, gpu_sigmaSq, gpu_n, gpu_epsilon_Cn);
}

__device__ double CalcVirShiftGPU(double distSq, int index, double gpu_sigmaSq,
                                  double *gpu_n, double *gpu_epsilon_Cn) {
  double rNeg2 = 1.0 / distSq;
  double rRat2 = gpu_sigmaSq * rNeg2;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] * 0.5);
  // return gpu_epsilon_Cn[index] * 6.0 *
  // ((gpu_n[index] / 6.0) * repulse - attract) * rNeg2;
  return gpu_epsilon_Cn[index] * rNeg2 *
         (gpu_n[index] * repulse - 6.0 * attract);
}

__device__ double CalcVirExp6GPU(double distSq, int index, double gpu_sigmaSq,
                                 double *gpu_n, double *gpu_rMin,
                                 double *gpu_rMaxSq, double *gpu_expConst,
                                 double sc_sigma_6, double sc_alpha,
                                 uint sc_power, double gpu_lambdaVDW) {
  if (distSq < gpu_rMaxSq[index]) {
    return DBL_MAX;
  }
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcVirExp6GPU(distSq, index, gpu_n, gpu_rMin, gpu_expConst);
  }

  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  return gpu_lambdaVDW * correction * correction *
         CalcVirExp6GPU(softRsq, index, gpu_n, gpu_rMin, gpu_expConst);
}

__device__ double CalcVirExp6GPU(double distSq, int index, double *gpu_n,
                                 double *gpu_rMin, double *gpu_expConst) {
  double dist = sqrt(distSq);
  double rRat = gpu_rMin[index] / dist;
  double rRat2 = rRat * rRat;
  double attract = rRat2 * rRat2 * rRat2;

  uint alpha_ij = gpu_n[index];
  double repulse =
      (dist / gpu_rMin[index]) * exp(alpha_ij * (1.0 - dist / gpu_rMin[index]));
  return 6.0 * gpu_expConst[index] * (repulse - attract) / distSq;
}

__device__ double CalcVirSwitchMartiniGPU(double distSq, int index,
                                          double gpu_sigmaSq, double *gpu_n,
                                          double *gpu_epsilon_Cn,
                                          double gpu_rCut, double gpu_rOn,
                                          double sc_sigma_6, double sc_alpha,
                                          uint sc_power, double gpu_lambdaVDW) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcVirSwitchMartiniGPU(distSq, index, gpu_sigmaSq, gpu_n,
                                   gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
  }

  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  return gpu_lambdaVDW * correction * correction *
         CalcVirSwitchMartiniGPU(softRsq, index, gpu_sigmaSq, gpu_n,
                                 gpu_epsilon_Cn, gpu_rCut, gpu_rOn);
}

__device__ double CalcVirSwitchMartiniGPU(double distSq, int index,
                                          double gpu_sigmaSq, double *gpu_n,
                                          double *gpu_epsilon_Cn,
                                          double gpu_rCut, double gpu_rOn) {
  double r_1 = rsqrt(distSq);
  double r_8 = distSq * distSq * distSq * distSq;
  double r_n2 = pow(r_1, gpu_n[index] + 2.0);

  double rij_ron = sqrt(distSq) - gpu_rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;

  double pn = gpu_n[index];
  double An =
      pn * ((pn + 1.0) * gpu_rOn - (pn + 4.0) * gpu_rCut) /
      (pow(gpu_rCut, pn + 2.0) * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double Bn = -pn * ((pn + 1.0) * gpu_rOn - (pn + 3.0) * gpu_rCut) /
              (pow(gpu_rCut, pn + 2.0) * (gpu_rCut - gpu_rOn) *
               (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));

  double sig6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  double sign = pow(gpu_sigmaSq, pn * 0.5);

  double A6 =
      6.0 * (7.0 * gpu_rOn - 10.0 * gpu_rCut) /
      (pow(gpu_rCut, 8.0) * (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));
  double B6 = -6.0 * (7.0 * gpu_rOn - 9.0 * gpu_rCut) /
              (pow(gpu_rCut, 8.0) * (gpu_rCut - gpu_rOn) *
               (gpu_rCut - gpu_rOn) * (gpu_rCut - gpu_rOn));

  double dshifttempRep = An * rij_ron_2 + Bn * rij_ron_3;
  double dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;

  const double dshiftRep =
      (distSq > gpu_rOn * gpu_rOn ? dshifttempRep * r_1 : 0);
  const double dshiftAtt =
      (distSq > gpu_rOn * gpu_rOn ? dshifttempAtt * r_1 : 0);
  double Wij = gpu_epsilon_Cn[index] * (sign * (pn * r_n2 + dshiftRep) -
                                        sig6 * (6.0 * r_8 + dshiftAtt));
  return Wij;
}

__device__ double CalcVirSwitchGPU(double distSq, int index, double gpu_sigmaSq,
                                   double *gpu_epsilon_Cn, double *gpu_n,
                                   double gpu_rCut, double gpu_rOn,
                                   double sc_sigma_6, double sc_alpha,
                                   uint sc_power, double gpu_lambdaVDW) {
  if (gpu_lambdaVDW >= 0.999999) {
    return CalcVirSwitchGPU(distSq, index, gpu_sigmaSq, gpu_epsilon_Cn, gpu_n,
                            gpu_rCut, gpu_rOn);
  }

  double sigma6 = gpu_sigmaSq * gpu_sigmaSq * gpu_sigmaSq;
  sigma6 = max(sigma6, sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = sc_alpha * pow((1.0 - gpu_lambdaVDW), (double)sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  return gpu_lambdaVDW * correction * correction *
         CalcVirSwitchGPU(softRsq, index, gpu_sigmaSq, gpu_epsilon_Cn, gpu_n,
                          gpu_rCut, gpu_rOn);
}

__device__ double CalcVirSwitchGPU(double distSq, int index, double gpu_sigmaSq,
                                   double *gpu_epsilon_Cn, double *gpu_n,
                                   double gpu_rCut, double gpu_rOn) {
  double rCutSq = gpu_rCut * gpu_rCut;
  double rCutSq_rijSq = rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;
  double rOnSq = gpu_rOn * gpu_rOn;

  double rNeg2 = 1.0 / distSq;
  double rRat2 = rNeg2 * gpu_sigmaSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double repulse = pow(rRat2, gpu_n[index] * 0.5);
  double factor1 = rCutSq - 3.0 * rOnSq;
  double factor2 =
      1.0 / ((rCutSq - rOnSq) * (rCutSq - rOnSq) * (rCutSq - rOnSq));

  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2.0 * distSq);
  double fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

  const double factE = (distSq > rOnSq ? fE : 1.0);
  const double factW = (distSq > rOnSq ? fW : 0.0);

  double Wij =
      gpu_epsilon_Cn[index] * rNeg2 * (gpu_n[index] * repulse - 6.0 * attract);
  double Eij = gpu_epsilon_Cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

#endif
