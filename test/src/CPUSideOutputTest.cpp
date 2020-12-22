#include <gtest/gtest.h>
#include "CPUSide.h"
#include "Simulation.h"

TEST(CPUSideOutput, CheckForConsistencyPDB_PSF_DCD) {

    std::string inputFileString = "outputTest.conf";
    Simulation sim(inputFileString.c_str());
    
    CPUSide * cpu;

    

}
