if(ENSEMBLE_NVT)
   add_executable(NVT ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(NVT PROPERTIES 
      OUTPUT_NAME ${NVT_name}
      COMPILE_FLAGS "${NVT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NVT ws2_32)
   endif()
   if(MPI_FOUND)
      target_link_libraries(NVT ${MPI_LIBRARIES})
      target_link_libraries(NVT stdc++fs)
   endif()
endif()

if(ENSEMBLE_GEMC)
   add_executable(GEMC ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(GEMC PROPERTIES 
      OUTPUT_NAME ${GE_name}
      COMPILE_FLAGS "${GE_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GEMC ws2_32)
   endif()
   if(MPI_FOUND)
      target_link_libraries(GEMC ${MPI_LIBRARIES})
      target_link_libraries(GEMC stdc++fs)
   endif()
endif()

if(ENSEMBLE_GCMC)
   add_executable(GCMC ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(GCMC PROPERTIES 
      OUTPUT_NAME ${GC_name}
      COMPILE_FLAGS "${GC_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GCMC ws2_32)
   endif()
   if(MPI_FOUND)
      target_link_libraries(GCMC ${MPI_LIBRARIES})
      target_link_libraries(GCMC stdc++fs)
   endif()
endif()

if(ENSEMBLE_NPT)
   add_executable(NPT ${sources} ${headers} ${libHeaders} ${libSources})
   set_target_properties(NPT PROPERTIES 
      OUTPUT_NAME ${NPT_name}
      COMPILE_FLAGS "${NPT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NPT ws2_32)
   endif()
   if(MPI_FOUND)
      target_link_libraries(NPT ${MPI_LIBRARIES})
      target_link_libraries(NPT stdc++fs)
   endif()
endif()

