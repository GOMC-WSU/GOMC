cmake_minimum_required(VERSION 3.18)

project(GOMC)
# set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(lib)
include_directories(src)
include_directories(src/cbmc)
include_directories(src/moves)
include_directories(src/GPU)

# Versioning
set(GOMC_VERSION_MAJOR 2)
set(GOMC_VERSION_MINOR 75)

# Out-of-source build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Run clang-tidy
if(GOMC_TIDY)
   set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,mpi-*,openmp-*)
endif()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Set compile and link flags. Need to do it this way so that we can pass
# the source compiler flags to NVCC properly
set(CMAKE_INTEL_COMP_FLAGS -O3 -fp-model=precise)
set(CMAKE_INTEL_CUDA_COMP_FLAGS "SHELL:-Xcompiler -O3" "SHELL:-fp-model=precise")
set(CMAKE_INTEL_LINK_FLAGS -O3 -fp-model=precise)
set(CMAKE_GNU_COMP_FLAGS -O3 -m64)
set(CMAKE_GNU_CUDA_COMP_FLAGS "SHELL:-Xcompiler -O3" "SHELL:-Xcompiler -m64")
set(CMAKE_GNU_LINK_FLAGS -m64)
set(CMAKE_CLANG_COMP_FLAGS -O3)
set(CMAKE_CLANG_CUDA_COMP_FLAGS "SHELL:-Xcompiler -O3")
set(CMAKE_CLANG_LINK_FLAGS -O3)

# The GOMC_OPT setting enables processor-specific optimizations
if(GOMC_OPT)
   set(CMAKE_INTEL_COMP_FLAGS ${CMAKE_INTEL_COMP_FLAGS} -xHost)
   set(CMAKE_INTEL_CUDA_COMP_FLAGS ${CMAKE_INTEL_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -xHost")
   set(CMAKE_INTEL_LINK_FLAGS ${CMAKE_INTEL_LINK_FLAGS} -xHost)
   set(CMAKE_GNU_COMP_FLAGS ${CMAKE_GNU_COMP_FLAGS} -march=native)
   set(CMAKE_GNU_CUDA_COMP_FLAGS ${CMAKE_GNU_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -march=native")
   set(CMAKE_GNU_LINK_FLAGS ${CMAKE_GNU_LINK_FLAGS} -march=native)
   set(CMAKE_CLANG_COMP_FLAGS ${CMAKE_CLANG_COMP_FLAGS} -march=native)
   set(CMAKE_CLANG_CUDA_COMP_FLAGS ${CMAKE_CLANG_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -march=native")
   set(CMAKE_CLANG_LINK_FLAGS ${CMAKE_CLANG_LINK_FLAGS} -march=native)
endif()


# Set the necessary flags for using address sanitizer, which doesn't work with the Intel compiler
if(GOMC_ASAN)
   set(CMAKE_GNU_COMP_FLAGS ${CMAKE_GNU_COMP_FLAGS} -fsanitize=address -fno-omit-frame-pointer)
   set(CMAKE_GNU_CUDA_COMP_FLAGS ${CMAKE_GNU_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -fsanitize=address" "SHELL:-Xcompiler -fno-omit-frame-pointer")
   set(CMAKE_GNU_LINK_FLAGS ${CMAKE_GNU_LINK_FLAGS} -fsanitize=address -fno-omit-frame-pointer)
   set(CMAKE_CLANG_COMP_FLAGS ${CMAKE_CLANG_COMP_FLAGS} -fsanitize=address -fno-omit-frame-pointer)
   set(CMAKE_CLANG_CUDA_COMP_FLAGS ${CMAKE_CLANG_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -fsanitize=address" "SHELL:-Xcompiler -fno-omit-frame-pointer")
   set(CMAKE_CLANG_LINK_FLAGS ${CMAKE_CLANG_LINK_FLAGS} -fsanitize=address -fno-omit-frame-pointer)
   set(CMAKE_MVSC_ASAN_FLAGS "/fsanitize=address")
endif()

# Find OpenMP and set it up
# Don't use OpenMP if using address sanitizer
if(NOT GOMC_ASAN)
   find_package(OpenMP)
   if(OPENMP_FOUND)
      set(CMAKE_INTEL_COMP_FLAGS ${CMAKE_INTEL_COMP_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_INTEL_CUDA_COMP_FLAGS ${CMAKE_INTEL_CUDA_COMP_FLAGS} "SHELL:-Xcompiler ${OpenMP_CXX_FLAGS}")
      set(CMAKE_INTEL_LINK_FLAGS ${CMAKE_INTEL_LINK_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_GNU_COMP_FLAGS ${CMAKE_GNU_COMP_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_GNU_CUDA_COMP_FLAGS ${CMAKE_GNU_CUDA_COMP_FLAGS} "SHELL:-Xcompiler ${OpenMP_CXX_FLAGS}")
      set(CMAKE_GNU_LINK_FLAGS ${CMAKE_GNU_LINK_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_CLANG_COMP_FLAGS ${CMAKE_CLANG_COMP_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_CLANG_CUDA_COMP_FLAGS ${CMAKE_CLANG_CUDA_COMP_FLAGS} "SHELL:-Xcompiler ${OpenMP_CXX_FLAGS}")
      set(CMAKE_CLANG_LINK_FLAGS ${CMAKE_CLANG_LINK_FLAGS} ${OpenMP_CXX_FLAGS})
      set(CMAKE_MSVC_OPENMP_FLAGS "/openmp:llvm")
   endif()
endif()

set(ENSEMBLE_NVT ON CACHE BOOL "Build NVT version")
set(ENSEMBLE_GEMC ON CACHE BOOL "Build GEMC version")
set(ENSEMBLE_GCMC ON CACHE BOOL "Build GCMC version")
set(ENSEMBLE_NPT ON CACHE BOOL "Build NPT version")
set(ENSEMBLE_GPU_NVT ON CACHE BOOL "Build GPU NVT version")
set(ENSEMBLE_GPU_GEMC ON CACHE BOOL "Build GPU GEMC version")
set(ENSEMBLE_GPU_GCMC ON CACHE BOOL "Build GPU GCMC version")
set(ENSEMBLE_GPU_NPT ON CACHE BOOL "Build GPU NPT version")

include(${PROJECT_SOURCE_DIR}/CMake/GOMCMPI.cmake)

include_directories("${PROJECT_BINARY_DIR}")

# Set compiler flags now that all the options have been specified
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "IntelLLVM")
   set(CMAKE_COMP_FLAGS ${CMAKE_INTEL_COMP_FLAGS})
   set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_INTEL_CUDA_COMP_FLAGS} -allow-unsupported-compiler)
   set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler)
   set(CMAKE_LINK_FLAGS ${CMAKE_INTEL_LINK_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
   set(CMAKE_COMP_FLAGS ${CMAKE_GNU_COMP_FLAGS})
   set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_GNU_CUDA_COMP_FLAGS})
   set(CMAKE_LINK_FLAGS ${CMAKE_GNU_LINK_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
   set(CMAKE_COMP_FLAGS ${CMAKE_CLANG_COMP_FLAGS})
   set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CLANG_CUDA_COMP_FLAGS})
   set(CMAKE_LINK_FLAGS ${CMAKE_CLANG_LINK_FLAGS})
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
   set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_MSVC_OPENMP_FLAGS} ${CMAKE_MVSC_ASAN_FLAGS} /D DEBUG /MTd /Zi /Ob0 /Od /RTC1 /D_USE_MATH_DEFINES)
   set(CMAKE_CXX_FLAGS_DEBUG_INIT ${CMAKE_CXX_FLAGS_DEBUG_INIT} ${CMAKE_MSVC_OPENMP_FLAGS} ${CMAKE_MVSC_ASAN_FLAGS} /D DEBUG /MTd /Zi /Ob0 /Od /RTC1 /D_USE_MATH_DEFINES)
   set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT ${CMAKE_CXX_FLAGS_MINSIZEREL_INIT} ${CMAKE_MSVC_OPENMP_FLAGS} /MT /O1 /Ob1 /D NDEBUG /D_USE_MATH_DEFINES)
   set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_MINSIZEREL_INIT} ${CMAKE_MSVC_OPENMP_FLAGS} /MT /O2 /D NDEBUG /D_USE_MATH_DEFINES)
   set(CMAKE_CXX_FLAGS_RELEASE_INIT ${CMAKE_CXX_FLAGS_RELEASE_INIT} ${CMAKE_MSVC_OPENMP_FLAGS} /MT /O2 /D NDEBUG /D_USE_MATH_DEFINES)
   set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT ${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} ${CMAKE_MSVC_OPENMP_FLAGS} /MT /Zi /O2 /Ob1 /D NDEBUG /D_USE_MATH_DEFINES)
endif()

# Set Source and Header files
include(${PROJECT_SOURCE_DIR}/CMake/FileLists.cmake)

# Optional CMake Argument
option(GOMC_GTEST    "Build unit tests for GOMC" OFF)
option(GOMC_GTEST_MPI    "Build unit tests for GOMC - MPI Enabled" OFF)

# Enable config header
configure_file(
    "${PROJECT_SOURCE_DIR}/GOMC_Config.h.in"
    "${PROJECT_BINARY_DIR}/GOMC_Config.h"
)

# Enable google test
# for now we will disable testing for Intel compiler
if(GOMC_GTEST OR GOMC_GTEST_MPI)
  enable_testing()
  include(test/GoogleTest.cmake)
endif()

# Check if CUDA exists and, if so, set the host compiler and enable CUDA
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
   set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
   enable_language(CUDA)
   include(${PROJECT_SOURCE_DIR}/CMake/GOMCCUDASetup.cmake)
else()
   message(WARNING "No CUDA support for this compiler")
endif()

# Set up Serial version
include(${PROJECT_SOURCE_DIR}/CMake/GOMCCPUSetup.cmake)
