#include "CUDAMemoryManager.cuh"

#ifdef GOMC_CUDA
long long CUDAMemoryManager::totalAllocatedBytes = 0;
std::unordered_map<void *, std::pair<unsigned int, std::string> > CUDAMemoryManager::allocatedPointers = {};


cudaError_t CUDAMemoryManager::mallocMemory(void **address, unsigned int size, std::string var_name) {
  cudaError_t ret = cudaMalloc(address, size);
  allocatedPointers[*address] = make_pair(size, var_name);
  totalAllocatedBytes += size;
  return ret;
}

cudaError_t CUDAMemoryManager::freeMemory(void *address) {
  if(allocatedPointers.find(address) != allocatedPointers.end()) {
    totalAllocatedBytes -= allocatedPointers[address].first;
    allocatedPointers.erase(address);
  } else {
    std::cout << "Warning! You are trying to free memory where it was freed" <<
      "\tor never been allocated before!\n";
  }
  return cudaFree(address);
}

bool CUDAMemoryManager::isFreed() {
  bool ret = allocatedPointers.size() == 0;
  while(allocatedPointers.size() != 0) {
    auto it = allocatedPointers.begin();
    std::cout << "You forgot to free memory address " << it->second.second
      << " with " << it->second.first << " bytes allocated to it!\n";
    std::cout << "I am going to free it for you!\n";
    freeMemory(it->first);
  }
  return ret;
}

#endif