#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <unordered_map>
#include <iostream>

class CUDAMemoryManager {
public:
  static cudaError_t mallocMemory(void **address, unsigned int size);
  static cudaError_t freeMemory(void *address);

private:
  static long long totalAllocatedBytes;
  static std::unordered_map<void *, unsigned int> allocatedPointers;
  static bool isFreed();
};

cudaError_t CUDAMemoryManager::mallocMemory(void **address, unsigned int size) {
  cudaError_t ret = cudaMalloc(address, size);
  allocatedPointers[*address] = size;
  totalAllocatedBytes += size;
  return ret;
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

bool CUDAMemoryManager::isFreed() {
  bool ret = allocatedPointers.size() == 0;
  while(allocatedPointers.size() != 0) {
    auto it = allocatedPointers.begin();
    std::cout << "You forgot to free memory address " << it->first
      << " with " << it->second << " bytes allocated to it!\n";
    std::cout << "I am going to free it for you!\n";
    freeMemory(it->first);
  }
  return ret;
}

#endif