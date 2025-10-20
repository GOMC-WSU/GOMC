# EnsemblePreprocessor defines:
# NVT = 1, GEMC = 2, GCMC = 3, NPT = 4

include(GoogleTest)

# Try to find MPI if GOMC_GTEST_MPI is ON
if(GOMC_GTEST_MPI)
  find_package(MPI REQUIRED)
endif()

function(add_ensemble_test name ENSEMBLE_ID ENSEMBLE_LABEL)
    set(flags "-DENSEMBLE=${ENSEMBLE_ID}")

    add_executable(${name}
        ${GOMCHeaders} ${GOMCSources}
        ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources}
    )

    # Link to gtest and optionally MPI
    target_link_libraries(${name} gtest_main)
    if(GOMC_GTEST_MPI)
      target_link_libraries(${name} MPI::MPI_CXX)
      target_compile_definitions(${name} PRIVATE USE_MPI)
    endif()

    set_target_properties(${name} PROPERTIES COMPILE_FLAGS "${flags}")

    add_test(NAME BasicTypesTest_${ENSEMBLE_LABEL}
             COMMAND ${name} --gtest_filter=BasicTypesTest.*)
    add_test(NAME MolLookupTest_${ENSEMBLE_LABEL}
             COMMAND ${name} --gtest_filter=CheckConsensusBeta.*)

    if(${ENSEMBLE_ID} EQUAL 2 OR ${ENSEMBLE_ID} EQUAL 3)
        add_test(NAME ConsistentTrajectoryTest_${ENSEMBLE_LABEL}
                 COMMAND ${name} --gtest_filter=CheckPDBTrajCoordinates.*)
    endif()
endfunction()

function(add_NVT_test name)
    add_ensemble_test(${name} 1 NVT)
endfunction()

function(add_NPT_test name)
    add_ensemble_test(${name} 4 NPT)
endfunction()

function(add_GCMC_test name)
    add_ensemble_test(${name} 3 GCMC)
endfunction()

function(add_GEMC_test name)
    add_ensemble_test(${name} 2 GEMC)
endfunction()

# Register all test executables
if(GOMC_GTEST_MPI)
  add_NVT_test(GOMC_NVT_MPI_Test)
  add_NPT_test(GOMC_NPT_MPI_Test)
  add_GCMC_test(GOMC_GCMC_MPI_Test)
  add_GEMC_test(GOMC_GEMC_MPI_Test)
else()
  add_NVT_test(GOMC_NVT_Test)
  add_NPT_test(GOMC_NPT_Test)
  add_GCMC_test(GOMC_GCMC_Test)
  add_GEMC_test(GOMC_GEMC_Test)
endif()
