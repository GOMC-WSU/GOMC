#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include "CalculateMinImageCUDAKernel.cuh"
#include<thrust/device_vector.h>
#include<thrust/sequence.h>
#include "CUDAMemoryManager.cuh"

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
    void MapParticlesToCell(VariablesCUDA * cv,
                            XYZArray const &coords,
                            XYZArray const &axes);
    void CopyGPUMemoryToToHost(int * deviceMemory,
                                    int size,
                                    std::vector<int> & hostMemory);
    void SortMappedParticles(VariablesCUDA * cv,
                              XYZArray const &coords);
    void CalculateCellDegrees(VariablesCUDA * cv,
                              XYZArray const &coords);
    void CalculateCellDegreesCUB(VariablesCUDA * cv,
                              XYZArray const &coords);
    void PrefixScanCellDegrees(VariablesCUDA * cv,
                                        int numberOfCells);
    // new pair interaction calculation done on GPU
    int *gpu_cellVector, *gpu_mapParticleToCell, *gpu_cellStartIndex;
    // Fixed as long as volume doesnt change
    // Regenerate after volume moves.
    int *gpu_neighborList, *gpu_numberOfCells, *gpu_startOfBoxCellList;
    int *gpu_edgeCells;
    double *gpu_cellSize;
  private:
    void CreateStartVector(int numberOfAtoms,
                          int * mapParticleToCell,
                          int * mapParticleToCellSorted,
                          int * particleIndices,
                          int * particleIndicesSorted);

    void CreateCellDegrees(int numberOfAtoms,
                          int * mapParticleToCellSortedGPURes,
                          int * OnesGPURes,
                          int * gpu_CellDegreeSanityCheck,
                          int * cellDegreesGPURes,
                          int * iterationsReq);

    void CalculateNewRowOffsets( int numberOfRows,
                                int * global_row_offsets_dev_ptr,
                                int * global_degrees_dev_ptr);
};

#endif
#endif