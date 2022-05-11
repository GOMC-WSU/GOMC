#ifdef GOMC_CUDA
#include "CellListGPU.cuh"

CellListGPU::CellListGPU(VariablesCUDA * cv):
gpu_cellVector(cv->gpu_cellVectorGPURes),
gpu_mapParticleToCell(cv->gpu_mapParticleToCellGPURes),
gpu_cellStartIndex(cv->gpu_cellStartIndexGPURes),
gpu_neighborList(cv->gpu_neighborListGPURes),
gpu_numberOfCells(cv->gpu_numberOfCellsGPURes),
gpu_startOfBoxCellList(cv->gpu_startOfBoxCellListGPURes)
{
}
#endif
