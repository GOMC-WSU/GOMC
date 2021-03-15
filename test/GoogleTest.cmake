function(add_mpi_test name no_mpi_proc)
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCMPISources} ${GOMCMPIHeaders} ${MPITestSources})
      add_dependencies(${name} googletest)
  # Make sure to link MPI here too:
  target_link_libraries(${name} ${MPI_LIBRARIES} gtest_main)
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

if(MPI_FOUND)
  add_mpi_test(GOMC_MPI_Test 2)
else()
# Now simply link against gtest or gtest_main as needed. Eg
add_executable(GOMC_Test ${sources} ${headers} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
target_link_libraries(GOMC_Test gtest_main)
add_test(NAME BasicTypesTest COMMAND BasicTypesTest)
add_test(NAME CircuitTester COMMAND DialaTest)
add_test(NAME MolLookupTest COMMAND CheckConsensusBeta)
add_test(NAME PSFParserTest COMMAND CheckProtAndWaterTest)
add_test(NAME EndianTest COMMAND TestBitSwap)
<<<<<<< 6999568989b504076e46b53df19131353683d601
add_test(NAME ConsistentTrajectoryTest COMMAND CheckPDBTrajCoordinates)


add_test(NAME ParallelTemperingTest COMMAND ParallelTemperingTest)
if(MPI_FOUND)
  target_link_libraries(GOMC_Test ${MPI_LIBRARIES})
endif()
=======
endif()




>>>>>>> Replica communication tested for com and coords
