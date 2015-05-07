cmake_minimum_required(VERSION 2.6)

project(GOMC_Serial)

#Out-of-source build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

#Versioning
set (GOMC_VERSION_MAJOR 0)
set (GOMC_VERSION_MINOR 97)

IF(NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(ENSEMBLE_NVT ON CACHE BOOL "Build NVT version")
set(ENSEMBLE_GEMC ON CACHE BOOL "Build GEMC version")
set(ENSEMBLE_GCMC ON CACHE BOOL "Build GCMC version")


#enable config header
configure_file(
	"${PROJECT_SOURCE_DIR}/GOMC_Config.h.in"
	"${PROJECT_BINARY_DIR}/GOMC_Config.h"
	)
include_directories("${PROJECT_BINARY_DIR}")

if(CMAKE_COMPILER_IS_CNUCXX)
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall") 
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall") 
endif()



#EnsemblePreprocessor defines NVT = 1, GEMC = 2, GCMC = 3
#Grand Canonical Monte Carlo
set(GC_flags "-DENSEMBLE=3")
set(GC_name "GOMC_Serial_GCMC")

#Gibbs Ensemble Monte Carlo
set(GE_flags "-DENSEMBLE=2")
set(GE_name "GOMC_Serial_GEMC")

#NVT (Canonical) Ensemble
set(NVT_flags "-DENSUMBLE=1")
set(NVT_name "GOMC_Serial_NVT")

add_subdirectory(src)
