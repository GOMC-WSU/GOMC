#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"
#include "XYZArray.h"
#include "CalculateMinImageCUDAKernel.cuh"

class CellListGPU {
  public:
    CellListGPU(VariablesCUDA * cv);
    void MapParticlesToCell(VariablesCUDA * cv,
                            XYZArray const &coords,
                            XYZArray const &axes);
    // new pair interaction calculation done on GPU
    int *gpu_cellVector, *gpu_mapParticleToCell, *gpu_cellStartIndex;
    // Fixed as long as volume doesnt change
    // Regenerate after volume moves.
    int *gpu_neighborList, *gpu_numberOfCells, *gpu_startOfBoxCellList;
    int *gpu_cellSize, *gpu_edgeCells;
  private:
    __device__ int PositionToCell(int atomIndex,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            double* gpu_cellSize,
                            int* gpu_edgeCells,
                            int* gpu_nonOrth,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z);

};

#endif
#endif