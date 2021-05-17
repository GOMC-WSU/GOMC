cmake_minimum_required(VERSION 3.8)

project(GOMC)

include_directories(lib)
include_directories(src)
include_directories(src/cbmc)
include_directories(src/moves)
include_directories(src/GPU)

# Enable google test
# for now we will disable testing for intel compiler
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  enable_testing()
  include(test/GoogleTest.cmake)
endif()

#Out-of-source build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#Versioning
set (GOMC_VERSION_MAJOR 2)
set (GOMC_VERSION_MINOR 70)

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

include(${PROJECT_SOURCE_DIR}/CMake/GOMCMPI.cmake)

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

# Set Source and Header files
include(${PROJECT_SOURCE_DIR}/CMake/FileLists.cmake)

# find Boost, and if not found download and install serialization library for checkpointing
# Set Source and Header files
#set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/FindBoost.cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
find_package(Boost COMPONENTS serialization )

if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
endif()

include_directories(${BOOST_INCLUDE_DIRS})

# Find if CUDA exists and what is the version number
include(CheckLanguage)
check_language(CUDA)
if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    include(${PROJECT_SOURCE_DIR}/CMake/GOMCCUDASetup.cmake)
endif()

# Setup Serial version
include(${PROJECT_SOURCE_DIR}/CMake/GOMCCPUSetup.cmake)

# find OpenMP and set it up
find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


