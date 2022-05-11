#ifdef GOMC_CUDA
#include "CellListGPU.cuh"
#include "CalculateMinImageCUDAKernel.cuh"

CellListGPU::CellListGPU(VariablesCUDA * cv):
gpu_cellVector(cv->gpu_cellVectorGPURes),
gpu_mapParticleToCell(cv->gpu_mapParticleToCellGPURes),
gpu_cellStartIndex(cv->gpu_cellStartIndexGPURes),
gpu_neighborList(cv->gpu_neighborListGPURes),
gpu_numberOfCells(cv->gpu_numberOfCellsGPURes),
gpu_startOfBoxCellList(cv->gpu_startOfBoxCellListGPURes),
gpu_cellSize(cv->gpu_cellSizeGPURes),
gpu_edgeCells(cv->gpu_edgeCellsGPURes)
{
}

void CellListGPU::MapParticlesToCell(VariablesCUDA * cv,
                                    XYZArray const &coords,
                                    XYZArray const &axes){
    int atomNumber = coords.Count();
    cudaMemcpy(cv->gpu_x, coords.x, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(cv->gpu_y, coords.y, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(cv->gpu_z, coords.z, atomNumber * sizeof(double), cudaMemcpyHostToDevice);
}

__device__ int PositionToCell(int atomIndex,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            double* gpu_cellSize,
                            int* gpu_edgeCells,
                            int* gpu_nonOrth,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z){
    double3 pos = make_double3(gpu_x[atomIndex],
                               gpu_y[atomIndex],
                               gpu_z[atomIndex]);
    if(gpu_nonOrth[0]){
        pos = TransformUnSlant(pos, 
                                gpu_Invcell_x,
                                gpu_Invcell_y,
                                gpu_Invcell_z);
    }
    int x = (int)(pos.x / gpu_cellSize[0]);
    int y = (int)(pos.y / gpu_cellSize[1]);
    int z = (int)(pos.z / gpu_cellSize[2]);
    //Check the cell number to avoid segfaults for coordinates close to axis
    //x, y, and z should never be equal or greater than number of cells in x, y,
    // and z axis, respectively.
    x -= (x == gpu_edgeCells[0] ?  1 : 0);
    y -= (y == gpu_edgeCells[1] ?  1 : 0);
    z -= (z == gpu_edgeCells[2] ?  1 : 0);
    return x * gpu_edgeCells[1] * gpu_edgeCells[2] + y * gpu_edgeCells[2] + z;
}

__global__ void MapParticlesToCell(int atomNumber,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            int* gpu_mapParticleToCell,
                            double *gpu_cellSize,
                            int *gpu_edgeCells,
                            int* gpu_nonOrth,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadID >= atomNumber)
        return;
    int cell = PositionToCell(threadID,
                            gpu_x,
                            gpu_y,
                            gpu_z,
                            gpu_cellSize,
                            gpu_edgeCells,
                            gpu_nonOrth,
                            gpu_Invcell_x,
                            gpu_Invcell_y,
                            gpu_Invcell_z);

}
                            

#endif
