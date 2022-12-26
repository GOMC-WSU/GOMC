# Find CUDA is enabled, set it up

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("-- Debug build type detected, passing : '-g -G --keep -lineinfo' to nvcc")
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -G --keep -lineinfo")
endif()


set(GEN_COMP_flag "-DGOMC_CUDA -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT ")

if (GOMC_NVTX_ENABLED)
	message("-- Enabling profiling with NVTX for GPU")
	set(GEN_COMP_flag "${GEN_COMP_flag} -DGOMC_NVTX_ENABLED")
endif()


include_directories(src/GPU)

set(GPU_NPT_flags "-DENSEMBLE=4 ${GEN_COMP_flag}")
set(GPU_NPT_name "GOMC_GPU_NPT")
set(GPU_GC_flags "-DENSEMBLE=3 ${GEN_COMP_flag}")
set(GPU_GC_name "GOMC_GPU_GCMC")
set(GPU_GE_flags "-DENSEMBLE=2 ${GEN_COMP_flag}")
set(GPU_GE_name "GOMC_GPU_GEMC")
set(GPU_NVT_flags "-DENSEMBLE=1 ${GEN_COMP_flag}")
set(GPU_NVT_name "GOMC_GPU_NVT")

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# Set host compiler
set(CCBIN "-ccbin=${CMAKE_CXX_COMPILER}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CCBIN} -Wno-deprecated-gpu-targets")

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

#####################################
function(add_GPU_NVT_test name)
    add_executable(${name} ${cudaSources} ${cudaHeaders} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources})
    target_link_libraries(${name} gtest_main ${GLOG_DEPENDENCY})
    set_target_properties(${name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${name}
        CUDA_ARCHITECTURES 35 60 70
        COMPILE_FLAGS "${GPU_NVT_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, ${name} setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET ${name} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(${name} ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(${name} ${MPI_LIBRARIES})
    endif()
endfunction(add_GPU_NVT_test)

function(add_GPU_NPT_test name)
    add_executable(${name} ${cudaSources} ${cudaHeaders} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources})
    target_link_libraries(${name} gtest_main ${GLOG_DEPENDENCY})
    set_target_properties(${name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${name}
        CUDA_ARCHITECTURES 35 60 70
        COMPILE_FLAGS "${GPU_NPT_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, ${name} setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET ${name} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(${name} ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(${name} ${MPI_LIBRARIES})
    endif()
endfunction(add_GPU_NPT_test)

function(add_GPU_GCMC_test name)
    add_executable(${name} ${cudaSources} ${cudaHeaders} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources})
    target_link_libraries(${name} gtest_main ${GLOG_DEPENDENCY})
    set_target_properties(${name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${name}
        CUDA_ARCHITECTURES 35 60 70
        COMPILE_FLAGS "${GPU_GC_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, ${name} setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET ${name} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(${name} ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(${name} ${MPI_LIBRARIES})
    endif()
endfunction(add_GPU_GCMC_test)

function(add_GPU_GEMC_test name)
    add_executable(${name} ${cudaSources} ${cudaHeaders} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources})
    target_link_libraries(${name} gtest_main ${GLOG_DEPENDENCY})
    set_target_properties(${name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${name}
        CUDA_ARCHITECTURES 35 60 70
        COMPILE_FLAGS "${GPU_GE_flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, ${name} setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET ${name} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    if(WIN32)
        target_link_libraries(${name} ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(${name} ${MPI_LIBRARIES})
    endif()
endfunction(add_GPU_GEMC_test)

add_GPU_NVT_test(GOMC_GPU_NVT_Test)
add_GPU_NPT_test(GOMC_GPU_NPT_Test)
add_GPU_GCMC_test(GOMC_GPU_GCMC_Test)
add_GPU_GEMC_test(GOMC_GPU_GEMC_Test)
set(GOMC_GTEST 1)
