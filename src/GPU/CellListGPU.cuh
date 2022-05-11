#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"

class CellListGPU {
  public:
    CellListGPU(VariablesCUDA * cv);
    void GridAll(VariablesCUDA * cv);
    // new pair interaction calculation done on GPU
    int *gpu_cellVector, *gpu_mapParticleToCell, *gpu_cellStartIndex;
    // Fixed as long as volume doesnt change
    // Regenerate after volume moves.
    int *gpu_neighborList, *gpu_numberOfCells, *gpu_startOfBoxCellList;
  private:
    void MapParticlesToCell(int atomNumber,
                            int* gpu_mapParticleToCell);
};

#endif
#endif