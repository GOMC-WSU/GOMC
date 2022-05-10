#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"

class CellListGPU {
    CellListGPU(VariablesCUDA * cv);
  // new pair interaction calculation done on GPU
  int *gpu_cellVector;
  int *gpu_mapParticleToCell;

};

#endif
#endif