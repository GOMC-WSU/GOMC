# HIP/ROCm build setup for GOMC

# Define flags for both CXX and HIP compilation of the GPU targets only, so
# the CPU ensembles in the same build stay CPU-only.
# __HIP_PLATFORM_AMD__ is needed when CXX files include HIP headers
set(GPU_HIP_defs GOMC_CUDA USE_HIP __HIP_PLATFORM_AMD__)

# enable_language(HIP) in the top-level CMakeLists already resolved
# CMAKE_HIP_ARCHITECTURES: it honors an explicit -DCMAKE_HIP_ARCHITECTURES,
# otherwise auto-detects the host GPU(s) and errors on a no-GPU host. Collapse
# to unique values before it is embedded in the custom --offload-arch link rule
# below (auto-detection lists an arch once per GPU agent).
list(REMOVE_DUPLICATES CMAKE_HIP_ARCHITECTURES)

# On Windows with Clang, add WIN32 (source guards) and _USE_MATH_DEFINES (M_2_SQRTPI etc.)
# Also override the HIP device-link rule: CMake 4.x Windows-Clang injects
# -fuse-ld=lld-link into <LINK_FLAGS>, which the --hip-link device-link mode
# rejects as an invalid linker name. Route through hip_link_win.py, which
# strips -fuse-ld=* and converts MSVC-style linker flags for the GCC driver.
if(WIN32)
    add_compile_definitions(WIN32 _USE_MATH_DEFINES)
    set(_HIP_LINK_WRAP "${CMAKE_CURRENT_SOURCE_DIR}/CMake/hip_link_win.py")
    set(CMAKE_HIP_LINK_EXECUTABLE
        "python3 ${_HIP_LINK_WRAP} ${CMAKE_HIP_COMPILER} --offload-arch=${CMAKE_HIP_ARCHITECTURES} <OBJECTS> -fgpu-rdc --hip-link -o <TARGET> <LINK_FLAGS> <LINK_LIBRARIES>")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("-- Debug build type detected for HIP")
endif()

message(STATUS "HIP architectures: ${CMAKE_HIP_ARCHITECTURES}")

include_directories(src/GPU)

set(GPU_NVT_flags "-DENSEMBLE=1")
set(GPU_NVT_name "GOMC_GPU_NVT")
set(GPU_GE_flags "-DENSEMBLE=2")
set(GPU_GE_name "GOMC_GPU_GEMC")
set(GPU_GC_flags "-DENSEMBLE=3")
set(GPU_GC_name "GOMC_GPU_GCMC")
set(GPU_NPT_flags "-DENSEMBLE=4")
set(GPU_NPT_name "GOMC_GPU_NPT")

set(CMAKE_HIP_STANDARD 17)
set(CMAKE_HIP_STANDARD_REQUIRED true)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# Mark CUDA sources as HIP language and enable RDC
set_source_files_properties(${cudaSources} PROPERTIES LANGUAGE HIP)
# Add RDC flags for device code linking across translation units
# Use default visibility to allow device functions to be linked across TUs
set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -fgpu-rdc")

# Find hipCUB - use include dirs directly to avoid flag pollution from hip::device
find_package(hipcub REQUIRED)
# Get include dirs without propagating compile flags
get_target_property(HIPCUB_INCLUDE_DIRS hip::hipcub INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${HIPCUB_INCLUDE_DIRS})
# Find rocprim for includes
find_package(rocprim REQUIRED)
get_target_property(ROCPRIM_INCLUDE_DIRS roc::rocprim INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${ROCPRIM_INCLUDE_DIRS})

#####################################
if(ENSEMBLE_GPU_NVT)
    add_executable(GPU_NVT ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    target_compile_options(GPU_NVT
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:HIP>:${CMAKE_COMP_FLAGS}>)
    target_link_options(GPU_NVT
       PUBLIC ${CMAKE_LINK_FLAGS})
    set_target_properties(GPU_NVT PROPERTIES
        HIP_SEPARABLE_COMPILATION ON
        HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}"
        OUTPUT_NAME ${GPU_NVT_name})
    target_compile_definitions(GPU_NVT PRIVATE ${GPU_HIP_defs} ENSEMBLE=1)
    target_link_libraries(GPU_NVT PRIVATE amdhip64)
    if(WIN32)
        target_link_libraries(GPU_NVT PRIVATE ws2_32)
    endif()
    if(MPI_FOUND)
        target_link_libraries(GPU_NVT ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GEMC)
    add_executable(GPU_GEMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    target_compile_options(GPU_GEMC
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:HIP>:${CMAKE_COMP_FLAGS}>)
    target_link_options(GPU_GEMC
       PUBLIC ${CMAKE_LINK_FLAGS})
    set_target_properties(GPU_GEMC PROPERTIES
        HIP_SEPARABLE_COMPILATION ON
        HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}"
        OUTPUT_NAME ${GPU_GE_name})
    target_compile_definitions(GPU_GEMC PRIVATE ${GPU_HIP_defs} ENSEMBLE=2)
    target_link_libraries(GPU_GEMC PRIVATE amdhip64)
    if(WIN32)
        target_link_libraries(GPU_GEMC PRIVATE ws2_32)
    endif()
    if(MPI_FOUND)
        target_link_libraries(GPU_GEMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_GCMC)
    add_executable(GPU_GCMC ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    target_compile_options(GPU_GCMC
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:HIP>:${CMAKE_COMP_FLAGS}>)
    target_link_options(GPU_GCMC
       PUBLIC ${CMAKE_LINK_FLAGS})
    set_target_properties(GPU_GCMC PROPERTIES
        HIP_SEPARABLE_COMPILATION ON
        HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}"
        OUTPUT_NAME ${GPU_GC_name})
    target_compile_definitions(GPU_GCMC PRIVATE ${GPU_HIP_defs} ENSEMBLE=3)
    target_link_libraries(GPU_GCMC PRIVATE amdhip64)
    if(WIN32)
        target_link_libraries(GPU_GCMC PRIVATE ws2_32)
    endif()
    if(MPI_FOUND)
        target_link_libraries(GPU_GCMC ${MPI_LIBRARIES})
    endif()
endif()

if(ENSEMBLE_GPU_NPT)
    add_executable(GPU_NPT ${cudaSources} ${cudaHeaders} ${sources}
    ${headers} ${libHeaders} ${libSources})
    target_compile_options(GPU_NPT
       PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${CMAKE_COMP_FLAGS}>
              $<$<COMPILE_LANGUAGE:HIP>:${CMAKE_COMP_FLAGS}>)
    target_link_options(GPU_NPT
       PUBLIC ${CMAKE_LINK_FLAGS})
    set_target_properties(GPU_NPT PROPERTIES
        HIP_SEPARABLE_COMPILATION ON
        HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}"
        OUTPUT_NAME ${GPU_NPT_name})
    target_compile_definitions(GPU_NPT PRIVATE ${GPU_HIP_defs} ENSEMBLE=4)
    target_link_libraries(GPU_NPT PRIVATE amdhip64)
    if(WIN32)
        target_link_libraries(GPU_NPT PRIVATE ws2_32)
    endif()
    if(MPI_FOUND)
        target_link_libraries(GPU_NPT ${MPI_LIBRARIES})
    endif()
endif()
