#EnsemblePreprocessor defines NVT = 1, GEMC = 2, GCMC = 3, NPT = 4
#NPT (Isothermal-Isobaric) Ensemble

set(NPT_flags "-DENSEMBLE=4")
set(NPT_name "GOMC_CPU_NPT")

#Grand Canonical Monte Carlo
set(GC_flags "-DENSEMBLE=3")
set(GC_name "GOMC_CPU_GCMC")

#Gibbs Ensemble Monte Carlo
set(GE_flags "-DENSEMBLE=2")
set(GE_name "GOMC_CPU_GEMC")

#NVT (Canonical) Ensemble
set(NVT_flags "-DENSEMBLE=1")
set(NVT_name "GOMC_CPU_NVT")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)

set(LIBRARIES_TO_LINK "${GLOG_DEPENDENCY}")
if(WIN32)
   set(LIBRARIES_TO_LINK "${LIBRARIES_TO_LINK} ws2_32")
endif()
if(MPI_FOUND)
   set(LIBRARIES_TO_LINK "${LIBRARIES_TO_LINK} ${MPI_LIBRARIES}")
endif()

if(ENSEMBLE_NVT)
   add_executable(NVT ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(NVT PROPERTIES 
      OUTPUT_NAME ${NVT_name}
      COMPILE_FLAGS "${NVT_flags}")
   target_link_libraries(NVT ${LIBRARIES_TO_LINK})
endif()

if(ENSEMBLE_GEMC)
   add_executable(GEMC ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(GEMC PROPERTIES 
      OUTPUT_NAME ${GE_name}
      COMPILE_FLAGS "${GE_flags}")
   target_link_libraries(GEMC ${LIBRARIES_TO_LINK})
endif()

if(ENSEMBLE_GCMC)
   add_executable(GCMC ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(GCMC PROPERTIES 
      OUTPUT_NAME ${GC_name}
      COMPILE_FLAGS "${GC_flags}")
   target_link_libraries(GCMC ${LIBRARIES_TO_LINK})
endif()

if(ENSEMBLE_NPT)
   add_executable(NPT ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(NPT PROPERTIES 
      OUTPUT_NAME ${NPT_name}
      COMPILE_FLAGS "${NPT_flags}")
   target_link_libraries(NPT ${LIBRARIES_TO_LINK})
endif()

