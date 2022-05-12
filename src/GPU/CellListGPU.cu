#ifdef GOMC_CUDA
#include "CellListGPU.cuh"
#include "cub/cub.cuh"



CellListGPU::CellListGPU(VariablesCUDA * cv, int atomNumber)
{
    CUMALLOC((void**) &cv->gpu_Ones, atomNumber * sizeof(int));
    CUMALLOC((void**) &cv->gpu_particleIndices, atomNumber * sizeof(int));
    thrust::device_vector<int>pI(atomNumber);
	thrust::sequence(pI.begin(),pI.end());
    thrust::device_vector<int>ones(atomNumber);
    thrust::fill(ones.begin(),ones.end(), 1);
	thrust::sequence(pI.begin(),pI.end());
    // Fill with 1s
    cudaMemcpy(cv->gpu_Ones, thrust::raw_pointer_cast(&ones[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
    cudaMemcpy(cv->gpu_particleIndices, thrust::raw_pointer_cast(&pI[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
}

void CellListGPU::GridAll(VariablesCUDA * cv,
                        XYZArray const &coords,
                        XYZArray const &axes,
                        int numberOfCells){
    GOMC_EVENT_START(1, GomcProfileEvent::GRID_ALL_GPU);
    MapParticlesToCell(cv,coords,axes);
    SortMappedParticles(cv,coords);
    CalculateCellDegrees(cv,coords);
    PrefixScanCellDegrees(cv, numberOfCells);
    GOMC_EVENT_STOP(1, GomcProfileEvent::GRID_ALL_GPU);
}


void CellListGPU::MapParticlesToCell(VariablesCUDA * cv,
                                    XYZArray const &coords,
                                    XYZArray const &axes){
    int atomNumber = coords.Count();
    // Run the kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (int)(atomNumber / threadsPerBlock) + 1;
    MapParticlesToCellKernel<<< blocksPerGrid, threadsPerBlock>>>(
                            atomNumber,
                            cv->gpu_x,
                            cv->gpu_y,
                            cv->gpu_z,                                
                            cv->gpu_mapParticleToCell,
                            cv->gpu_cellSize,
                            cv->gpu_edgeCells,
                            cv->gpu_nonOrth,
                            cv->gpu_Invcell_x[0],
                            cv->gpu_Invcell_y[0],
                            cv->gpu_Invcell_z[0]);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}


void CellListGPU::SortMappedParticles(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();
    // Run the kernel
    CreateStartVector(atomNumber,
                    cv->gpu_mapParticleToCell,
                    cv->gpu_mapParticleToCellSorted,
                    cv->gpu_particleIndices,
                    cv->gpu_cellVector,
                    cv->d_temp_storage_sort,
                    cv->temp_storage_bytes_sort);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}

void CellListGPU::CalculateCellDegrees(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();                                            // Run the kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (int)(atomNumber / threadsPerBlock) + 1;
    // Run the kernel
    /*
    CalculateCellDegreesKernel<<< blocksPerGrid, 
                                threadsPerBlock, 
                                (3*threadsPerBlock+2)*sizeof(int)>>>(atomNumber,
                                                                cv->gpu_mapParticleToCellSorted,
                                                                cv->gpu_cellDegrees);
    */
 CalculateCellDegreesKernel<<< blocksPerGrid, 
                                threadsPerBlock>>>(atomNumber,
                                                    cv->gpu_mapParticleToCellSorted,
                                                    cv->gpu_cellDegrees);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}

void CellListGPU::CalculateCellDegreesCUB(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();                                            // Run the kernel
    CreateCellDegrees(atomNumber,
                    cv->gpu_mapParticleToCellSorted,
                    cv->gpu_Ones,
                    cv->gpu_CellDegreeSanityCheck,
                    cv->gpu_cellDegrees,
                    cv->gpu_IterationsReq);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}
// CustomMin functor
struct CustomMin
{
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__
    T operator()(const T &a, const T &b) const {
        return (b < a) ? b : a;
    }
};
// Need to sort first, since ReduceByKey needs keys to be in stretches of the same key
void CellListGPU::CreateCellDegrees(int numberOfAtoms,
                                int * mapParticleToCellSorted,
                                int * Ones,
                                int * gpu_CellDegreeSanityCheck,
                                int * cellDegrees,
                                int * iterationsReq){
    // Declare, allocate, and initialize device-accessible pointers for input and output
    int          num_items = numberOfAtoms;          // e.g., 8
    int          *d_keys_in = mapParticleToCellSorted;           // e.g., [0, 2, 2, 9, 5, 5, 5, 8]
    int          *d_values_in = Ones;                            // e.g., [0, 7, 1, 6, 2, 5, 3, 4]
    int          *d_unique_out = gpu_CellDegreeSanityCheck;      // e.g., [-, -, -, -, -, -, -, -]
    int          *d_aggregates_out = cellDegrees;  // e.g., [-, -, -, -, -, -, -, -]
    int          *d_num_runs_out = iterationsReq;    // e.g., [-]
    CustomMin    reduction_op;
    // Determine temporary device storage requirements
    void     *d_temp_storage = NULL;
    size_t   temp_storage_bytes = 0;
    cub::DeviceReduce::ReduceByKey(d_temp_storage, temp_storage_bytes, d_keys_in, d_unique_out, d_values_in, d_aggregates_out, d_num_runs_out, reduction_op, num_items);
    // Allocate temporary storage
    cudaMalloc(&d_temp_storage, temp_storage_bytes);
    // Run reduce-by-key
    cub::DeviceReduce::ReduceByKey(d_temp_storage, temp_storage_bytes, d_keys_in, d_unique_out, d_values_in, d_aggregates_out, d_num_runs_out, reduction_op, num_items);
    // d_unique_out      <-- [0, 2, 9, 5, 8]
    // d_aggregates_out  <-- [0, 1, 6, 2, 4]
    // d_num_runs_out    <-- [5]
}


void CellListGPU::PrefixScanCellDegrees(VariablesCUDA * cv,
                                    int numberOfCells){
    CalculateNewRowOffsets(numberOfCells,
                           cv->gpu_cellStartIndex,
                           cv->gpu_cellDegrees,
                           cv->d_temp_storage_prefix_sum,
                           cv->temp_storage_bytes_prefix_sum);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);
}

void CellListGPU::CopyGPUMemoryToToHost(int * deviceMemory,
                                    int size,
                                    std::vector<int> & hostMemory){
    hostMemory.clear();
    hostMemory.resize(size);
    cudaMemcpy(&hostMemory[0], deviceMemory, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);
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
        pos = TransformUnSlantGPU(pos, 
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
    gpu_mapParticleToCell[threadID] = cell;

}



// Make sure a zero element is padded onto the end.
// https://github.com/NVIDIA/cub/issues/367
void CellListGPU::CalculateNewRowOffsets( int numberOfRows,
                                        int * global_row_offsets_dev_ptr,
                                        int * global_degrees_dev_ptr,
                                        void     **d_temp_storage_prefix_sum,
                                        size_t   *temp_storage_bytes_prefix_sum){
    // Declare, allocate, and initialize device-accessible pointers for input and output
    int  num_items = numberOfRows+1;      // e.g., 7
    int  *d_in = global_degrees_dev_ptr;        // e.g., [8, 6, 7, 5, 3, 0, 9]
    int  *d_out = global_row_offsets_dev_ptr;         // e.g., [ ,  ,  ,  ,  ,  ,  ]
    // Determine temporary device storage requirements
    if(temp_storage_bytes_prefix_sum == NULL){
        cub::DeviceScan::ExclusiveSum(*d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum, d_in, d_out, num_items);
        // Allocate temporary storage
        CUMALLOC(d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum);
    } else {
        cudaMemset(*d_temp_storage_prefix_sum, 0, *temp_storage_bytes_prefix_sum);
    }
    // Run exclusive prefix sum
    cub::DeviceScan::ExclusiveSum(*d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum, d_in, d_out, num_items);
    // d_out s<-- [0, 8, 14, 21, 26, 29, 29]
    //cudaFree(temp_storage_bytes_prefix_sum);
}



void CellListGPU::CreateStartVector(int numberOfAtoms,
                                int * mapParticleToCell,
                                int * mapParticleToCellSorted,
                                int * particleIndices,
                                int * particleIndicesSorted,
                                void     **d_temp_storage_sort,
                                size_t   *temp_storage_bytes_sort){
    // Declare, allocate, and initialize device-accessible pointers for sorting data
    // numberOfAtoms            e.g., 7
    // mapParticleToCell        e.g., [8, 6, 7, 5, 3, 0, 9]
    // mapParticleToCellSorted  e.g., [        ...        ]
    // particleIndices          e.g., [0, 1, 2, 3, 4, 5, 6]
    // particleIndicesSorted    e.g., [        ...        ]
    // Determine temporary device storage requirements
    int num_items = numberOfAtoms;
    int  *d_keys_in = mapParticleToCell;
    int  *d_keys_out = mapParticleToCellSorted;       
    int  *d_values_in = particleIndices;    
    int  *d_values_out = particleIndicesSorted;   
    if(temp_storage_bytes_sort == NULL){
        cub::DeviceRadixSort::SortPairs(*d_temp_storage_sort, *temp_storage_bytes_sort,
            d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);
        // Allocate temporary storage
        CUMALLOC(d_temp_storage_sort, *temp_storage_bytes_sort);
    } else {
        cudaMemset(*d_temp_storage_sort, 0, *temp_storage_bytes_sort);
    }
    // Run sorting operation
    cub::DeviceRadixSort::SortPairs(*d_temp_storage_sort, *temp_storage_bytes_sort,
        d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);
    // mapParticleToCellSorted        <-- [0, 3, 5, 6, 7, 8, 9]
    // particleIndicesSorted          <-- [5, 4, 3, 1, 2, 0, 6]
}

__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadID >= atomNumber)
        return;
    atomicAdd(&gpu_cellDegrees[gpu_mapParticleToCellSorted[threadID]], 1);
}

/*
__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadID >= atomNumber)
        return;
    extern __shared__ int part_ary[];
    // minimum cell to maximum cell
    // 0..BlockSize-1 - count of each cell key
    // BlockSize..2*BlockSize-1 - indicator variable (0,1)
    // 2*BlockSize..3*BlockSize-1 - load cell key from global memory
    // 3*BlockSize - Min cell key
    // 3*BlockSize+1 - Max cell key
    part_ary[2*blockDim.x+threadIdx.x] = gpu_mapParticleToCellSorted[threadID];
    part_ary[threadIdx.x] = part_ary[2*blockDim.x+threadIdx.x];
    int i = blockDim.x/2;
    while (0 < i){
        if(threadIdx.x < i){
            part_ary[threadIdx.x] = min(part_ary[threadIdx.x],part_ary[threadIdx.x+i]); 
        }
          __syncthreads();
        i /= 2;
    }
    if (threadIdx.x==0)
        part_ary[3*blockDim.x] = part_ary[threadIdx.x];

      __syncthreads();

    part_ary[threadIdx.x] = part_ary[2*blockDim.x+threadIdx.x];
    i = blockDim.x/2;
    while (0 < i){
        if(threadIdx.x < i){
            part_ary[threadIdx.x] = max(part_ary[threadIdx.x],part_ary[threadIdx.x+i]); 
        }
          __syncthreads();
        i /= 2;
    }
    if (threadIdx.x==0)
        part_ary[3*blockDim.x+1] = part_ary[threadIdx.x];

      __syncthreads();

    for (int reductionIteration = 0; reductionIteration < part_ary[3*blockDim.x+1]-part_ary[3*blockDim.x]; ++reductionIteration){
        part_ary[blockDim.x+threadIdx.x] =  part_ary[2*blockDim.x+threadIdx.x] == part_ary[3*blockDim.x]+reductionIteration;
        i = blockDim.x/2;
        while (0 < i){
            if(threadIdx.x < i){
                part_ary[blockDim.x+threadIdx.x] = part_ary[blockDim.x+threadIdx.x] + part_ary[blockDim.x+threadIdx.x+i]; 
            }
              __syncthreads();
            i /= 2;
        }
        if (threadIdx.x==0)
            part_ary[reductionIteration] = part_ary[blockDim.x+threadIdx.x];

          __syncthreads();
    }

    for (int reductionIteration = threadIdx.x; reductionIteration < part_ary[3*blockDim.x+1]-part_ary[3*blockDim.x]; ++reductionIteration){
        atomicAdd(&gpu_cellDegrees[part_ary[3*blockDim.x]+reductionIteration], part_ary[reductionIteration]);
    }
}
*/

#endif
