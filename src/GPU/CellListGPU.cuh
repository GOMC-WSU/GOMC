#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include "CalculateMinImageCUDAKernel.cuh"
#include<thrust/device_vector.h>
#include<thrust/sequence.h>
#include<thrust/fill.h>
#include "CUDAMemoryManager.cuh"
#include "GOMCEventsProfile.h" // for NVTX profiling

// For cuMemsetD32
//#include <cuda_runtime.h>

__global__ void MapParticlesToCellKernel(int atomNumber,
                    double* gpu_x,
                    double* gpu_y,
                    double* gpu_z,                                
                    int* gpu_mapParticleToCell,
                    double *gpu_cellSize,
                    int *gpu_edgeCells,
                    int* gpu_nonOrth,
                    double *gpu_Invcell_x,
                    double *gpu_Invcell_y,
                    double *gpu_Invcell_z);

__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees);

class CellListGPU {
  public:
    CellListGPU(VariablesCUDA * cv, int atomCount);
    void GridAll(VariablesCUDA * cv,
                  XYZArray const &coords,
                  XYZArray const &axes,
                  int numberOfCells);
    void CopyGPUMemoryToToHost(int * deviceMemory,
                                    int size,
                                    std::vector<int> & hostMemory);
    void MapParticlesToCell(VariablesCUDA * cv,
                            XYZArray const &coords,
                            XYZArray const &axes);
    void SortMappedParticles(VariablesCUDA * cv,
                              XYZArray const &coords);
    void CalculateCellDegrees(VariablesCUDA * cv,
                              XYZArray const &coords);
    void CalculateCellDegreesCUB(VariablesCUDA * cv,
                              XYZArray const &coords);
    void PrefixScanCellDegrees(VariablesCUDA * cv,
                                        int numberOfCells);
  private:
    void CreateStartVector(int numberOfAtoms,
                          int * mapParticleToCell,
                          int * mapParticleToCellSorted,
                          int * particleIndices,
                          int * particleIndicesSorted,
                          void     *d_temp_storage_sort,
                          size_t   *temp_storage_bytes_sort);

    void CreateCellDegrees(int numberOfAtoms,
                          int * mapParticleToCellSortedGPURes,
                          int * OnesGPURes,
                          int * gpu_CellDegreeSanityCheck,
                          int * cellDegreesGPURes,
                          int * iterationsReq);

    void CalculateNewRowOffsets( int numberOfRows,
                                int * global_row_offsets_dev_ptr,
                                int * global_degrees_dev_ptr,
                                void     *d_temp_storage_prefix_sum,
                                size_t   *temp_storage_bytes_prefix_sum);
};

#endif
#endif