/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA

#include "CUDAMemoryManager.cuh"

long long CUDAMemoryManager::totalAllocatedBytes = 0;
std::unordered_map<void *, std::pair<unsigned int, std::string>>
    CUDAMemoryManager::allocatedPointers;

cudaError_t CUDAMemoryManager::mallocMemory(void **address, unsigned int size,
                                            std::string var_name) {
  cudaError_t ret = cudaMalloc(address, size);
  allocatedPointers[*address] = make_pair(size, var_name);
  totalAllocatedBytes += size;
  if (size == 0) {
    std::cout << "Warning! You are trying to allocate " << var_name
              << " with a size of zero bytes!\n";
  }
  return ret;
}

cudaError_t CUDAMemoryManager::freeMemory(void *address, std::string var_name) {
  if (allocatedPointers.find(address) != allocatedPointers.end()) {
    totalAllocatedBytes -= allocatedPointers[address].first;
    allocatedPointers.erase(address);
  } else if (address != nullptr) {
    std::cout << "Warning! You are trying to free " << var_name
              << " but it has already been freed\n"
              << "\tor was never allocated!\n";
  }
  return cudaFree(address);
}

bool CUDAMemoryManager::isFreed() {
  bool ret = allocatedPointers.size() == 0;
  while (allocatedPointers.size() != 0) {
    auto it = allocatedPointers.begin();
    std::cout << "You forgot to free memory " << it->second.second << " with "
              << it->second.first << " bytes allocated!\n";
    std::cout << "I am going to free it for you!\n";
    freeMemory(it->first, it->second.second);
  }
  return ret;
}

#endif
