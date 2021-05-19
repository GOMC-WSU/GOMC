set(TestSources
    test/src/BasicTypesTest.cpp
    test/src/BitLibTest.cpp
    test/src/EndianTest.cpp
    test/src/CircuitTester.cpp
    test/src/PSFParserTest.cpp
    test/src/ConsistentTrajectoryTest.cpp
)

set(TestHeaders

)

set(GOMCSources
    lib/AlphaNum.cpp
    src/PDBSetup.cpp
    src/MolSetup.cpp
    src/BondAdjacencyList.cpp
    src/ConfigSetup.cpp
    src/FFSetup.cpp
    src/FFConst.cpp
    src/Reader.cpp
    lib/CircuitFinder.cpp
    src/InputFileReader.cpp
    lib/FloydWarshallCycle.cpp
)

set(GOMCHeaders
    src/PDBSetup.h
    src/MolSetup.h
    lib/AlphaNum.h
    src/BondAdjacencyList.h
    src/ConfigSetup.h
    src/FFSetup.h
    src/FFConst.h
    src/Reader.h
    lib/CircuitFinder.h
    src/InputFileReader.h
    lib/FloydWarshallCycle.h
)