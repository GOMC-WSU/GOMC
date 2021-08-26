#EnsemblePreprocessor defines NVT = 1, GEMC = 2, GCMC = 3, NPT = 4

function(add_NVT_mpi_test name no_mpi_proc)
  set(NVT_flags "-DENSEMBLE=1")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCMPISources} ${GOMCMPIHeaders} ${MPITestSources})
      add_dependencies(${name} googletest)
  # Make sure to link MPI here too:
  target_link_libraries(${name} ${MPI_LIBRARIES} gtest_main)
  #set(test_parameters " ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc} ./${name}")
  #add_test(NAME ${name} COMMAND "${MPIEXEC} ${test_parameters}")
  add_test(NAME ${name} COMMAND ParallelTemperingTest)
endfunction(add_NVT_mpi_test)

function(add_NPT_mpi_test name no_mpi_proc)
  #NPT (Isothermal-Isobaric) Ensemble
  set(NPT_flags "-DENSEMBLE=4")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCMPISources} ${GOMCMPIHeaders} ${libHeaders} ${libSources} ${MPITestSources})
  # Make sure to link MPI here too:
  target_link_libraries(${name} ${MPI_LIBRARIES} gtest_main)
  set_target_properties(${name} PROPERTIES 
      COMPILE_FLAGS "${NPT_flags}")
  #set(test_parameters " ${MPIEXEC_NUMPROC_FLAG} ${no_mpi_proc} ./${name}")
  #add_test(NAME ${name} COMMAND "${MPIEXEC} ${test_parameters}")
  add_test(NAME ${name} COMMAND ParallelTemperingTest)
endfunction(add_mpi_test)

# Download and unpack googletest at configure time
configure_file(${PROJECT_SOURCE_DIR}/test/CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Include file lists
include(test/FileList.cmake)

if(GOMC_GTEST_MPI)
  add_NVT_mpi_test(GOMC_NVT_MPI_Test 2)
  add_NPT_mpi_test(GOMC_NPT_MPI_Test 2)
  set(GOMC_GTEST_MPI 1)
endif()

if(GOMC_GTEST)
# Now simply link against gtest or gtest_main as needed. Eg
add_executable(GOMC_Test ${sources} ${headers} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
target_link_libraries(GOMC_Test gtest_main)
add_test(NAME BasicTypesTest COMMAND BasicTypesTest)
add_test(NAME CircuitTester COMMAND DialaTest)
add_test(NAME MolLookupTest COMMAND CheckConsensusBeta)
add_test(NAME PSFParserTest COMMAND CheckProtAndWaterTest)
add_test(NAME EndianTest COMMAND TestBitSwap)
add_test(NAME ConsistentTrajectoryTest COMMAND CheckPDBTrajCoordinates)
endif()
