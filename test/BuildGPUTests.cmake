# Find CUDA is enabled, set it up
set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} -DGOMC_CUDA)
set(CMAKE_HOST_COMP_FLAGS ${CMAKE_COMP_FLAGS} -DGOMC_CUDA)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Debug build type detected, passing '-g -lineinfo --keep' to nvcc")
    set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -rdynamic" -g -lineinfo --keep)
    set(CMAKE_CUDA_LINK_FLAGS "SHELL:-Xcompiler -rdynamic" -g -lineinfo --keep)
endif()

if(GOMC_NVTX_ENABLED)
    message("-- Enabling profiling with NVTX for GPU")
    set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} -lineinfo -DGOMC_NVTX_ENABLED)
    set(CMAKE_HOST_COMP_FLAGS ${CMAKE_HOST_COMP_FLAGS} -DGOMC_NVTX_ENABLED)
endif()

# Set architecture flags based on the CMake version
# Once CMake 3.23 has been available for a while, we should just use
# set(CMAKE_CUDA_ARCHITECTURES all) and remove the if block
if (CMAKE_MAJOR_VERSION VERSION_GREATER 3 OR CMAKE_MINOR_VERSION VERSION_GREATER_EQUAL 23)
    set(CMAKE_CUDA_ARCHITECTURES all)
else()
# set the CUDA architectures based on the version of CUDA being used
    execute_process(COMMAND python ../scripts/get_cuda_full_version.py OUTPUT_VARIABLE NVCC_VERSION)
	if(NVCC_VERSION LESS 11.0)
        set(CMAKE_CUDA_ARCHITECTURES 35;37;50;52;53;60;61;70)
	elseif (NVCC_VERSION LESS 11.1)
        set(CMAKE_CUDA_ARCHITECTURES 52;53;60;61;70;75;80)
	elseif(NVCC_VERSION LESS 11.8)
        set(CMAKE_CUDA_ARCHITECTURES 52;53;60;61;70;75;80;86)
	elseif(NVCC_VERSION LESS 12.8)
        set(CMAKE_CUDA_ARCHITECTURES 52;53;60;61;70;75;80;86;89;90)
	else()
        set(CMAKE_CUDA_ARCHITECTURES 75;80;86;89;90;100;120)
	endif()
endif()

set(GEN_COMP_flag "-DGOMC_CUDA -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT ")

if (GOMC_NVTX_ENABLED)
	message("-- Enabling profiling with NVTX for GPU")
	set(GEN_COMP_flag "${GEN_COMP_flag} -DGOMC_NVTX_ENABLED")
endif()


include_directories(src/GPU)

set(GPU_NVT_flags "-DENSEMBLE=1")
set(GPU_GE_flags "-DENSEMBLE=2")
set(GPU_GC_flags "-DENSEMBLE=3")
set(GPU_NPT_flags "-DENSEMBLE=4")

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

#####################################
#####################################
function(add_GPU_ensemble_test name ENSEMBLE_ID)
    set(flags "-DENSEMBLE=${ENSEMBLE_ID}")
    add_executable(${name} ${cudaSources} ${cudaHeaders} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources})
    # Set compiler and linker flags for NVCC and the host compiler
    target_compile_options(${name}
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_HOST_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CUDA_COMP_FLAGS}>)
    target_link_options(${name}
       PUBLIC $<HOST_LINK:${CMAKE_LINK_FLAGS}>
              $<DEVICE_LINK:${CMAKE_CUDA_LINK_FLAGS}>)
    set_target_properties(${name} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${name}
        COMPILE_FLAGS "${flags}")
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		message("-- Debug build type detected, ${name} setting CUDA_RESOLVE_DEVICE_SYMBOLS ON")
    	set_property(TARGET ${name} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	endif()
    target_link_libraries(${name} gtest_main ${FFTW3_LIBRARIES})
    target_include_directories(${name} PRIVATE ${FFTW3_INCLUDE_DIRS})
    if(WIN32)
        target_link_libraries(${name} ws2_32)
    endif()
    if(MPI_FOUND)
	    target_link_libraries(${name} ${MPI_LIBRARIES})
    endif()
endfunction(add_GPU_ensemble_test)

function(add_GPU_NVT_test name)
    add_GPU_ensemble_test(${name} 1)
endfunction()

function(add_GPU_NPT_test name)
    add_GPU_ensemble_test(${name} 4)
endfunction()

function(add_GPU_GCMC_test name)
    add_GPU_ensemble_test(${name} 3)
endfunction()

function(add_GPU_GEMC_test name)
    add_GPU_ensemble_test(${name} 2)
endfunction()

# Register all test executables
if(GOMC_GTEST_MPI)
  add_GPU_NVT_test(GOMC_GPU_NVT_MPI_Test)
  add_GPU_NPT_test(GOMC_GPU_NPT_MPI_Test)
  add_GPU_GCMC_test(GOMC_GPU_GCMC_MPI_Test)
  add_GPU_GEMC_test(GOMC_GPU_GEMC_MPI_Test)
else()
  add_GPU_NVT_test(GOMC_GPU_NVT_Test)
  add_GPU_NPT_test(GOMC_GPU_NPT_Test)
  add_GPU_GCMC_test(GOMC_GPU_GCMC_Test)
  add_GPU_GEMC_test(GOMC_GPU_GEMC_Test)
endif()

set(GOMC_GTEST 1)
