# EnsemblePreprocessor defines:
# NVT = 1, GEMC = 2, GCMC = 3, NPT = 4

include(GoogleTest)  # Recommended, even if not using gtest_discover_tests()

# Helper function to reduce duplication
function(add_ensemble_test name ENSEMBLE_ID ENSEMBLE_LABEL)
    set(flags "-DENSEMBLE=${ENSEMBLE_ID}")

    # My tests are all called name_test.cpp
    add_executable(${name} 
        ${GOMCHeaders} ${GOMCSources}
        ${libHeaders} ${libSources}
        ${TestHeaders} ${TestSources}
    )

    target_link_libraries(${name} gtest_main)
    set_target_properties(${name} PROPERTIES COMPILE_FLAGS "${flags}")

    #
    # Register each suite as a separate CTest test
    # Each gtest_filter pattern limits execution to one test suite
    #
    add_test(NAME BasicTypesTest_${ENSEMBLE_LABEL}
             COMMAND ${name} --gtest_filter=BasicTypesTest.*)

    # add_test(NAME CircuitTester_${ENSEMBLE_LABEL}
    #          COMMAND ${name} --gtest_filter=DialaTest.*)

    add_test(NAME MolLookupTest_${ENSEMBLE_LABEL}
             COMMAND ${name} --gtest_filter=CheckConsensusBeta.*)

    # add_test(NAME PSFParserTest_${ENSEMBLE_LABEL}
    #          COMMAND ${name} --gtest_filter=CheckProtAndWaterTest.*)

    if(${ENSEMBLE_ID} EQUAL 2 OR ${ENSEMBLE_ID} EQUAL 3)  # GEMC or GCMC
        add_test(NAME ConsistentTrajectoryTest_${ENSEMBLE_LABEL}
                 COMMAND ${name} --gtest_filter=CheckPDBTrajCoordinates.*)
        # add_test(NAME CheckpointTest_${ENSEMBLE_LABEL}
        #          COMMAND ${name} --gtest_filter=CheckMollookup.*)
    endif()
endfunction()

# Define per-ensemble wrappers
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
add_NVT_test(GOMC_NVT_Test)
add_NPT_test(GOMC_NPT_Test)
add_GCMC_test(GOMC_GCMC_Test)
add_GEMC_test(GOMC_GEMC_Test)
