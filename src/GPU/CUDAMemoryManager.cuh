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
  static bool isFreed();

private:
  static long long totalAllocatedBytes;
  static std::unordered_map<void *, unsigned int> allocatedPointers;
};

#endif