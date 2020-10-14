########################  MPI   ##################################

include(${PROJECT_SOURCE_DIR}/CMake/GOMCOptionUtilities.cmake)

option(GOMC_MPI    "Build a parallel (message-passing) version of GOMC" OFF)
option(GOMC_THREAD_MPI  "Build a thread-MPI-based multithreaded version of GOMC (not compatible with MPI)" ON)

include(${PROJECT_SOURCE_DIR}/CMake/GOMCManageMPI.cmake)
include(${PROJECT_SOURCE_DIR}/CMake/ThreadMPI.cmake)

set(GOMC_COMMON_LIBRARIES "")
GOMC_dependent_option(
    GOMC_MPI_IN_PLACE
    "Enable MPI_IN_PLACE for MPIs that have it defined"
    ON
    GOMC_MPI)
mark_as_advanced(GOMC_MPI_IN_PLACE)

########################  MPI   ##################################