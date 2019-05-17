cmake_minimum_required(VERSION 2.6)

project(GOMC)

include_directories(lib)
include_directories(src)
include_directories(src/cbmc)
include_directories(src/moves)

#fix new policy of cmake about FindOpenMP.cmake (will have to check later if they fixed this issue)
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

#Out-of-source build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#Versioning
set (GOMC_VERSION_MAJOR 2)
set (GOMC_VERSION_MINOR 40)

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

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -stdlib=libc++") 
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -stdlib=libc++")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}") 
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
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
set(NPT_name "GOMC_CPU_NPT")

#Grand Canonical Monte Carlo
set(GC_flags "-DENSEMBLE=3")
set(GC_name "GOMC_CPU_GCMC")

#Gibbs Ensemble Monte Carlo
set(GE_flags "-DENSEMBLE=2")
set(GE_name "GOMC_CPU_GEMC")

#NVT (Canonical) Ensemble
set(NVT_flags "-DENSEMBLE=1")
set(NVT_name "GOMC_CPU_NVT")

# Set Source and Header files
include(${PROJECT_SOURCE_DIR}/CMake/FileLists.cmake)

# Setup Serial version
include(${PROJECT_SOURCE_DIR}/CMake/GOMCCPUSetup.cmake)

# find CUDA and set it up
set(CUDA_SEPARABLE_COMPILATION ON)
find_package(CUDA)
IF(CUDA_FOUND)
    include(${PROJECT_SOURCE_DIR}/CMake/GOMCCUDASetup.cmake)
ENDIF(CUDA_FOUND)

# find OpenMP and set it up
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
