/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#pragma once

// CUDA-to-HIP compatibility header for GOMC.
// When building with HIP, this header aliases CUDA symbols to their HIP
// equivalents. On CUDA builds, it simply includes the CUDA runtime.

#if defined(USE_HIP) || defined(__HIP_PLATFORM_AMD__)

#include <hip/hip_runtime.h>

// Error handling
#define cudaError_t                     hipError_t
#define cudaSuccess                     hipSuccess
#define cudaGetLastError                hipGetLastError
#define cudaGetErrorString              hipGetErrorString

// Memory management
#define cudaMalloc                      hipMalloc
#define cudaFree                        hipFree
#define cudaMemcpy                      hipMemcpy
#define cudaMemset                      hipMemset
#define cudaMemcpyHostToDevice          hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost          hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice        hipMemcpyDeviceToDevice
#define cudaMemGetInfo                  hipMemGetInfo

// Synchronization
#define cudaDeviceSynchronize           hipDeviceSynchronize

// Device management
#define cudaGetDeviceCount              hipGetDeviceCount
#define cudaSetDevice                   hipSetDevice
#define cudaGetDeviceProperties         hipGetDeviceProperties
#define cudaDeviceProp                  hipDeviceProp_t
#define cudaDeviceGetAttribute          hipDeviceGetAttribute
#define cudaDevAttrMemoryClockRate      hipDeviceAttributeMemoryClockRate

// Profiler (disabled on HIP; profiling via rocprof instead)
#define cudaProfilerStart()             ((void)0)
#define cudaProfilerStop()              ((void)0)

// CUB debug macro (hipCUB doesn't have CubDebugExit, so provide a simple version)
#define CubDebugExit(e) do { \
    hipError_t __cub_err = (e); \
    if (__cub_err != hipSuccess) { \
        fprintf(stderr, "HIP error %d at %s:%d\n", __cub_err, __FILE__, __LINE__); \
        exit(__cub_err); \
    } \
} while(0)

#else  // CUDA path

#include <cuda.h>
#include <cuda_runtime.h>

#endif
