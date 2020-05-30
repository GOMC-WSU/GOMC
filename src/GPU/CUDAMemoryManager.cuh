#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <unordered_map>
#include <iostream>

class CUDAMemoryManager {
public:
  void Init() { totalAllocatedBytes = 0; }
  cudaError_t mallocMemory(void **address, unsigned int size);
  cudaError_t freeMemory(void *address);

private:
  static long long totalAllocatedBytes;
  static std::unordered_map<void *, unsigned int> allocatedPointers;
};

cudaError_t CUDAMemoryManager::mallocMemory(void **address, unsigned int size) {
  allocatedPointers[*address] = size;
  totalAllocatedBytes += size;
  return cudaMalloc(address, size);
}

cudaError_t CUDAMemoryManager::freeMemory(void *address) {
  if(allocatedPointers.find(address) != allocatedPointers.end()) {
    totalAllocatedBytes -= allocatedPointers[address];
    allocatedPointers.erase(address);
  } else {
    std::cout << "Warning! You are trying to free memory where it was freed or never been allocated before!\n";
  }
  return cudaFree(address);
}

#endif