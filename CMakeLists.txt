cmake_minimum_required(VERSION 2.6)

project(GOMC_Serial)

# CUDA 
find_package(CUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_30,code=sm_30)

include_directories(lib)
include_directories(src)
include_directories(src/cbmc)
include_directories(src/GPU)

#Out-of-source build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#Versioning
set (GOMC_VERSION_MAJOR 2)
set (GOMC_VERSION_MINOR 00)

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(ENSEMBLE_NVT ON CACHE BOOL "Build NVT version")
set(ENSEMBLE_GEMC ON CACHE BOOL "Build GEMC version")
set(ENSEMBLE_GCMC ON CACHE BOOL "Build GCMC version")
set(ENSEMBLE_NPT ON CACHE BOOL "Build NPT version")
set(ENSEMBLE_GPU_NVT ON CACHE BOOL "Build GPU NVT version")
set(ENSEMBLE_GPU_GEMC ON CACHE BOOL "Build GPU GEMC version")
set(ENSEMBLE_GPU_GCMC ON CACHE BOOL "Build GPU GCMC version")
set(ENSEMBLE_GPU_NPT ON CACHE BOOL "Build GPU NPT version")

#enable config header
configure_file(
	"${PROJECT_SOURCE_DIR}/GOMC_Config.h.in"
	"${PROJECT_BINARY_DIR}/GOMC_Config.h"
	)
include_directories("${PROJECT_BINARY_DIR}")

if(CMAKE_COMPILER_IS_CNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall") 
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall") 
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "${CMAKE_CXX_FLAGS_DEBUG_INIT} /D_DEBUG /MTd /Zi /Ob0 /Od /RTC1")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "${CMAKE_CXX_FLAGS_MINSIZEREL_INIT} /MT /O1 /Ob1 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_MINSIZEREL_INIT} /MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "${CMAKE_CXX_FLAGS_RELEASE_INIT} /MT /O2 /Ob2 /D NDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} /MT /Zi /O2 /Ob1 /D NDEBUG")
endif()

#EnsemblePreprocessor defines NVT = 1, GEMC = 2, GCMC = 3, NPT = 4

#NPT (Isothermal-Isobaric) Ensemble
set(NPT_flags "-DENSEMBLE=4")
set(NPT_name "GOMC_Serial_NPT")

set(GPU_NPT_flags "-DENSEMBLE=4 -DGOMC_CUDA")
set(GPU_NPT_name "GOMC_GPU_NPT")

#Grand Canonical Monte Carlo
set(GC_flags "-DENSEMBLE=3")
set(GC_name "GOMC_Serial_GCMC")

set(GPU_GC_flags "-DENSEMBLE=3 -DGOMC_CUDA")
set(GPU_GC_name "GOMC_GPU_GCMC")

#Gibbs Ensemble Monte Carlo
set(GE_flags "-DENSEMBLE=2")
set(GE_name "GOMC_Serial_GEMC")

set(GPU_GE_flags "DENSEMBLE=2 -DGOMC_CUDA")
set(GPU_GE_name "GOMC_GPU_GEMC")

#NVT (Canonical) Ensemble
set(NVT_flags "-DENSEMBLE=1")
set(NVT_name "GOMC_Serial_NVT")

set(GPU_NVT_flags "-DENSEMBLE=1 -DGOMC_CUDA")
set(GPU_NVT_name "GOMC_GPU_NVT")

add_subdirectory(src)
