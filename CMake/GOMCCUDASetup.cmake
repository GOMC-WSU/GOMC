# Find CUDA is enabled, set it up
set(CMAKE_CUDA_COMP_FLAGS -DGOMC_CUDA)
set(CMAKE_HOST_COMP_FLAGS ${CMAKE_COMP_FLAGS} -DGOMC_CUDA)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Debug build type detected, passing '-g -lineinfo --keep' to nvcc")
    set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} "SHELL:-Xcompiler -rdynamic" -g -lineinfo --keep)
    set(CMAKE_CUDA_LINK_FLAGS "SHELL:-Xcompiler -rdynamic" -g -lineinfo --keep)
endif()

if(GOMC_NVTX_ENABLED)
    message("-- Enabling profiling with NVTX for GPU")
    set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} -DGOMC_NVTX_ENABLED)
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

include_directories(src/GPU)

set(GPU_NVT_flags "-DENSEMBLE=1")
set(GPU_NVT_name "GOMC_GPU_NVT")
set(GPU_GE_flags "-DENSEMBLE=2")
set(GPU_GE_name "GOMC_GPU_GEMC")
set(GPU_GC_flags "-DENSEMBLE=3")
set(GPU_GC_name "GOMC_GPU_GCMC")
set(GPU_NPT_flags "-DENSEMBLE=4")
set(GPU_NPT_name "GOMC_GPU_NPT")

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# Disable the warning on deprecated GPU targets
set(CMAKE_CUDA_COMP_FLAGS ${CMAKE_CUDA_COMP_FLAGS} -Wno-deprecated-gpu-targets)
set(CMAKE_CUDA_LINK_FLAGS ${CMAKE_CUDA_LINK_FLAGS} -Wno-deprecated-gpu-targets)

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

#####################################
if(ENSEMBLE_GPU_NVT)
    add_executable(GPU_NVT ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    # Set compiler and linker flags for NVCC and the host compiler
    target_compile_options(GPU_NVT
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_HOST_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CUDA_COMP_FLAGS}>)
    target_link_options(GPU_NVT
       PUBLIC $<HOST_LINK:${CMAKE_LINK_FLAGS}>
              $<DEVICE_LINK:${CMAKE_CUDA_LINK_FLAGS}>)
    set_target_properties(GPU_NVT PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_NVT_name}
        COMPILE_FLAGS "${GPU_NVT_flags}")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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
    # Set compiler and linker flags for NVCC and the host compiler
    target_compile_options(GPU_GEMC
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_HOST_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CUDA_COMP_FLAGS}>)
    target_link_options(GPU_GEMC
       PUBLIC $<HOST_LINK:${CMAKE_LINK_FLAGS}>
              $<DEVICE_LINK:${CMAKE_CUDA_LINK_FLAGS}>)
    set_target_properties(GPU_GEMC PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_GE_name}
        COMPILE_FLAGS "${GPU_GE_flags}")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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
    # Set compiler and linker flags for NVCC and the host compiler
    target_compile_options(GPU_GCMC
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_HOST_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CUDA_COMP_FLAGS}>)
    target_link_options(GPU_GCMC
       PUBLIC $<HOST_LINK:${CMAKE_LINK_FLAGS}>
              $<DEVICE_LINK:${CMAKE_CUDA_LINK_FLAGS}>)
    set_target_properties(GPU_GCMC PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_GC_name}
        COMPILE_FLAGS "${GPU_GC_flags}")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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
    # Set compiler and linker flags for NVCC and the host compiler
    target_compile_options(GPU_NPT
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_HOST_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:CUDA>:${CMAKE_CUDA_COMP_FLAGS}>)
    target_link_options(GPU_NPT
       PUBLIC $<HOST_LINK:${CMAKE_LINK_FLAGS}>
              $<DEVICE_LINK:${CMAKE_CUDA_LINK_FLAGS}>)
    set_target_properties(GPU_NPT PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        OUTPUT_NAME ${GPU_NPT_name}
        COMPILE_FLAGS "${GPU_NPT_flags}")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
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
