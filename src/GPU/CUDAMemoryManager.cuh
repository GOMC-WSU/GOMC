/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CUDA_MEMORY_MANAGER
#define CUDA_MEMORY_MANAGER

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <unordered_map>

#define CUMALLOC(address, size)                                                \
  CUDAMemoryManager::mallocMemory(address, size, #address)
#define CUFREE(address) CUDAMemoryManager::freeMemory(address, #address)

class CUDAMemoryManager {
public:
  static cudaError_t mallocMemory(void **address, unsigned int size,
                                  std::string var_name);
  static cudaError_t freeMemory(void *address, std::string var_name);
  static bool isFreed();

private:
  static long long totalAllocatedBytes;
  static std::unordered_map<void *, std::pair<unsigned int, std::string>>
      allocatedPointers;
};

#endif /*GOMC_CUDA*/
#endif /*CUDA_MEMORY_MANAGER*/
