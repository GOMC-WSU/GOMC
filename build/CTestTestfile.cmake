# CMake generated Testfile for 
# Source directory: /home/greg/GOMC
# Build directory: /home/greg/GOMC/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test("BasicTypesTest" "BasicTypesTest")
set_tests_properties("BasicTypesTest" PROPERTIES  _BACKTRACE_TRIPLES "/home/greg/GOMC/test/GoogleTest.cmake;39;add_test;/home/greg/GOMC/test/GoogleTest.cmake;0;;/home/greg/GOMC/CMakeLists.txt;46;include;/home/greg/GOMC/CMakeLists.txt;0;")
add_test("CircuitTester" "DialaTest")
set_tests_properties("CircuitTester" PROPERTIES  _BACKTRACE_TRIPLES "/home/greg/GOMC/test/GoogleTest.cmake;40;add_test;/home/greg/GOMC/test/GoogleTest.cmake;0;;/home/greg/GOMC/CMakeLists.txt;46;include;/home/greg/GOMC/CMakeLists.txt;0;")
add_test("MolLookupTest" "CheckConsensusBeta")
set_tests_properties("MolLookupTest" PROPERTIES  _BACKTRACE_TRIPLES "/home/greg/GOMC/test/GoogleTest.cmake;41;add_test;/home/greg/GOMC/test/GoogleTest.cmake;0;;/home/greg/GOMC/CMakeLists.txt;46;include;/home/greg/GOMC/CMakeLists.txt;0;")
add_test("PSFParserTest" "CheckProtAndWaterTest")
set_tests_properties("PSFParserTest" PROPERTIES  _BACKTRACE_TRIPLES "/home/greg/GOMC/test/GoogleTest.cmake;42;add_test;/home/greg/GOMC/test/GoogleTest.cmake;0;;/home/greg/GOMC/CMakeLists.txt;46;include;/home/greg/GOMC/CMakeLists.txt;0;")
add_test("EndianTest" "TestBitSwap")
set_tests_properties("EndianTest" PROPERTIES  _BACKTRACE_TRIPLES "/home/greg/GOMC/test/GoogleTest.cmake;43;add_test;/home/greg/GOMC/test/GoogleTest.cmake;0;;/home/greg/GOMC/CMakeLists.txt;46;include;/home/greg/GOMC/CMakeLists.txt;0;")
subdirs("googletest-build")
