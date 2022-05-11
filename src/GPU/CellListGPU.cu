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

void CellListGPU::MapParticlesToCell(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();
    cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
}

__global__ void MapParticlesToCell(int atomNumber,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            int* gpu_mapParticleToCell){}

#endif
