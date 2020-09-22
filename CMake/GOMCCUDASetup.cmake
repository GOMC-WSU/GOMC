# Find CUDA is enabled, set it up

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("-- Debug build type detected, passing : '-g -G --keep' to nvcc")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -G --keep")
endif()

include_directories(src/GPU)

set(GPU_NPT_flags "-DENSEMBLE=4 -DGOMC_CUDA")
set(GPU_NPT_name "GOMC_GPU_NPT")
set(GPU_GC_flags "-DENSEMBLE=3 -DGOMC_CUDA")
set(GPU_GC_name "GOMC_GPU_GCMC")
set(GPU_GE_flags "-DENSEMBLE=2 -DGOMC_CUDA")
set(GPU_GE_name "GOMC_GPU_GEMC")
set(GPU_NVT_flags "-DENSEMBLE=1 -DGOMC_CUDA")
set(GPU_NVT_name "GOMC_GPU_NVT")

if(NOT DEFINED CMAKE_CUDA_STANDARD)
   set(CMAKE_CUDA_STANDARD 14)
   set(CMAKE_CUDA_STANDARD_REQUIRED true)
endif()

find_package(CUDALibs)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

#####################################
if(ENSEMBLE_GPU_NVT)
    add_executable(GPU_NVT ${cudaSources} ${cudaHeaders}
    ${sources} ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_NVT PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_NVT_name}
        COMPILE_FLAGS "${GPU_NVT_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, GPU_NVT setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET GPU_NVT PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(GPU_NVT ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(GPU_NVT ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GEMC)
    add_executable(GPU_GEMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_GEMC PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_GE_name}
        COMPILE_FLAGS "${GPU_GE_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, GPU_GEMC setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET GPU_GEMC PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(GPU_GEMC ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(GPU_GEMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GCMC)
    add_executable(GPU_GCMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_GCMC PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_GC_name}
        COMPILE_FLAGS "${GPU_GC_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, GPU_GCMC setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET GPU_GCMC PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(GPU_GCMC ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(GPU_GCMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_NPT)
    add_executable(GPU_NPT ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    set_target_properties(GPU_NPT PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_NPT_name}
        COMPILE_FLAGS "${GPU_NPT_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, GPU_NPT setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET GPU_NPT PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(GPU_NPT ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(GPU_NPT ${MPI_LIBRARIES})
    endif()
endif()
