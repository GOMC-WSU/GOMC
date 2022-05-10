#ifdef GOMC_CUDA
#include "CellListGPU.cuh"

CellListGPU::CellListGPU(VariablesCUDA * cv):
gpu_cellVector(cv->gpu_cellVector),
gpu_mapParticleToCell(cv->gpu_mapParticleToCell),
gpu_cellStartIndex(cv->gpu_cellStartIndex),
gpu_neighborList(cv->gpu_neighborList),
gpu_numberOfCells(cv->gpu_numberOfCells),
gpu_startOfBoxCellList(cv->gpu_startOfBoxCellList)
{
}
#endif
