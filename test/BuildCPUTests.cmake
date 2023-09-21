#EnsemblePreprocessor defines NVT = 1, GEMC = 2, GCMC = 3, NPT = 4

function(add_NVT_test name)
      set(NVT_flags "-DENSEMBLE=1")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
      target_link_libraries(${name} gtest_main)
      set_target_properties(${name} PROPERTIES 
      COMPILE_FLAGS "${NVT_flags}")
      add_test(NAME BasicTypesTest_NVT COMMAND BasicTypesTest)
      #add_test(NAME CircuitTester_NVT COMMAND DialaTest)
      add_test(NAME MolLookupTest_NVT COMMAND CheckConsensusBeta)
      #add_test(NAME PSFParserTest_NVT COMMAND CheckProtAndWaterTest)
endfunction(add_NVT_test)

function(add_NPT_test name)
      set(NPT_flags "-DENSEMBLE=4")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
      target_link_libraries(${name} gtest_main)
      set_target_properties(${name} PROPERTIES 
      COMPILE_FLAGS "${NPT_flags}")
      add_test(NAME BasicTypesTest_NPT COMMAND BasicTypesTest)
      #add_test(NAME CircuitTester_NPT COMMAND DialaTest)
      add_test(NAME MolLookupTest_NPT COMMAND CheckConsensusBeta)
      #add_test(NAME PSFParserTest_NPT COMMAND CheckProtAndWaterTest)
endfunction(add_NPT_test)

function(add_GCMC_test name)
      set(GCMC_flags "-DENSEMBLE=3")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
      target_link_libraries(${name} gtest_main)
      set_target_properties(${name} PROPERTIES 
      COMPILE_FLAGS "${GCMC_flags}")
      add_test(NAME BasicTypesTest_GCMC COMMAND BasicTypesTest)
      #add_test(NAME CircuitTester_GCMC COMMAND DialaTest)
      add_test(NAME MolLookupTest_GCMC COMMAND CheckConsensusBeta)
      #add_test(NAME PSFParserTest_GCMC COMMAND CheckProtAndWaterTest)
      add_test(NAME ConsistentTrajectoryTest_GCMC COMMAND CheckPDBTrajCoordinates)
      #add_test(NAME CheckpointTest_GCMC COMMAND CheckMollookup)
endfunction(add_GCMC_test)

function(add_GEMC_test name)
      set(GEMC_flags "-DENSEMBLE=2")
      # My test are all called name_test.cpp
      add_executable(${name} ${GOMCHeaders} ${GOMCSources} ${libHeaders} ${libSources} ${TestHeaders} ${TestSources})
      target_link_libraries(${name} gtest_main)
      set_target_properties(${name} PROPERTIES 
      COMPILE_FLAGS "${GEMC_flags}")
      add_test(NAME BasicTypesTest_GEMC COMMAND BasicTypesTest)
      #add_test(NAME CircuitTester_GEMC COMMAND DialaTest)
      add_test(NAME MolLookupTest_GEMC COMMAND CheckConsensusBeta)
      #add_test(NAME PSFParserTest_GEMC COMMAND CheckProtAndWaterTest)
      add_test(NAME ConsistentTrajectoryTest_GEMC COMMAND CheckPDBTrajCoordinates)
      #add_test(NAME CheckpointTest_GEMC COMMAND CheckMollookup)
endfunction(add_GEMC_test)

add_NVT_test(GOMC_NVT_Test)
add_NPT_test(GOMC_NPT_Test)
add_GCMC_test(GOMC_GCMC_Test)
add_GEMC_test(GOMC_GEMC_Test)