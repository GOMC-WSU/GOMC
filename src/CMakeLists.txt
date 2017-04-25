set(sources 
   BlockOutput.cpp
   BoxDimensions.cpp
   CBMC.cpp
   CellList.cpp
   ConfigSetup.cpp
   ConsoleOutput.cpp
   Coordinates.cpp
   CPUSide.cpp
   CalculateEnergy.cpp
   EnPartCntSampleOutput.cpp
   Ewald.cpp
   EwaldCached.cpp
   FFConst.cpp
   FFDihedrals.cpp
   FFParticle.cpp
   FFSetup.cpp
   Forcefield.cpp
   Geometry.cpp
   HistOutput.cpp
   InputFileReader.cpp
   Main.cc
   MoleculeKind.cpp
   MoleculeLookup.cpp
   Molecules.cpp
   MolSetup.cpp
   MoveConst.cpp
   MoveSettings.cpp
   NoEwald.cpp
   OutConst.cpp
   OutputVars.cpp
   PDBSetup.cpp
   PDBOutput.cpp
   PRNGSetup.cpp
   PSFOutput.cpp
   Reader.cpp
   Simulation.cpp
   StaticVals.cpp
   System.cpp
   cbmc/DCGraph.cpp
   cbmc/DCFreeHedron.cpp
   cbmc/DCLinkedHedron.cpp
   cbmc/DCHedron.cpp
   cbmc/DCLinear.cpp
   cbmc/DCLink.cpp
   cbmc/DCLinkNoDih.cpp
   cbmc/DCOnSphere.cpp
   cbmc/DCSingle.cpp
   cbmc/TrialMol.cpp)

set(headers
   BlockOutput.h
   BoxDimensions.h
   CalculateEnergy.h
   CBMC.h
   CellList.h
   Clock.h
   COM.h
   ConfigSetup.h
   ConsoleOutput.h
   ConstField.h
   Coordinates.h
   CoordinateSetup.h
   CPUSide.h
   EnergyTypes.h
   EnPartCntSampleOutput.h
   EnsemblePreprocessor.h
   Ewald.h
   EwaldCached.h  
   FFAngles.h
   FFBonds.h
   FFConst.h
   FFDihedrals.h
   FFParticle.h
   FFSetup.h
   FFShift.h
   FFSwitch.h
   FFSwitchMartini.h
   FixedWidthReader.h
   Forcefield.h
   FxdWidthWrtr.h
   Geometry.h
   HistOutput.h
   InputAbstracts.h
   InputFileReader.h
   IntraSwap.h
   MersenneTwister.h
   MoleculeKind.h
   MoleculeLookup.h
   MoleculeTransfer.h
   Molecules.h
   MolPick.h
   MolSetup.h
   MoveBase.h
   MoveConst.h
   MoveSettings.h
   NoEwald.h
   OutConst.h
   OutputAbstracts.h
   OutputVars.h
   PDBConst.h
   PDBOutput.h
   PDBSetup.h
   PRNG.h
   PRNGSetup.h
   PSFOutput.h
   Reader.h
   SeedReader.h
   Setup.h
   SimEventFrequency.h
   Simulation.h
   StaticVals.h
   SubdividedArray.h
   System.h
   TransformMatrix.h
   Writer.h
   XYZArray.h
   cbmc/DCComponent.h
   cbmc/DCData.h
   cbmc/DCFreeHedron.h
   cbmc/DCLinkedHedron.h
   cbmc/DCGraph.h
   cbmc/DCHedron.h
   cbmc/DCLinear.h
   cbmc/DCLink.h
   cbmc/DCLinkNoDih.h
   cbmc/DCOnSphere.h
   cbmc/DCSingle.h
   cbmc/TrialMol.h
   GPU/CalculateEnergyCUDAKernel.h)

set(libHeaders
   ../lib/BasicTypes.h
   ../lib/BitLib.h
   ../lib/GeomLib.h
   ../lib/NumLib.h
   ../lib/StrLib.h
   ../lib/StrStrmLib.h
   ../lib/VectorLib.h)

set(cudaSources
   GPU/CalculateEnergyCUDAKernel.cu)

source_group("Header Files" FILES ${headers})
source_group("Lib Headers" FILES ${libHeaders})

if(ENSEMBLE_NVT)
   add_executable(NVT ${sources} ${headers} ${libHeaders})
   set_target_properties(NVT PROPERTIES 
      OUTPUT_NAME ${NVT_name}
      COMPILE_FLAGS "${NVT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NVT ws2_32)
   endif()
endif()

if(ENSEMBLE_GEMC)
   add_executable(GEMC ${sources} ${headers} ${libHeaders})
   set_target_properties(GEMC PROPERTIES 
      OUTPUT_NAME ${GE_name}
      COMPILE_FLAGS "${GE_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GEMC ws2_32)
   endif()
endif()

if(ENSEMBLE_GCMC)
   add_executable(GCMC ${sources} ${headers} ${libHeaders})
   set_target_properties(GCMC PROPERTIES 
      OUTPUT_NAME ${GC_name}
      COMPILE_FLAGS "${GC_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GCMC ws2_32)
   endif()
endif()

if(ENSEMBLE_NPT)
   add_executable(NPT ${sources} ${headers} ${libHeaders})
   set_target_properties(NPT PROPERTIES 
      OUTPUT_NAME ${NPT_name}
      COMPILE_FLAGS "${NPT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NPT ws2_32)
   endif()
endif()

if(ENSEMBLE_GPU_NVT)
   add_executable(GPU_NVT ${sources} ${headers} ${libHeaders} ${cudaSources})
   set_target_properties(GPU_NVT PROPERTIES 
      OUTPUT_NAME ${GPU_NVT_name}
      COMPILE_FLAGS "${GPU_NVT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NVT ws2_32)
   endif()
endif()

if(ENSEMBLE_GPU_GEMC)
   add_executable(GPU_GEMC ${sources} ${headers} ${libHeaders} ${cudaSources})
   set_target_properties(GPU_GEMC PROPERTIES 
      OUTPUT_NAME ${GPU_GE_name}
      COMPILE_FLAGS "${GPU_GE_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GEMC ws2_32)
   endif()
endif()

if(ENSEMBLE_GPU_GCMC)
   cuda_add_executable(GPU_GCMC ${sources} ${headers} ${libHeaders} ${cudaSources})
   set_target_properties(GPU_GCMC PROPERTIES 
      OUTPUT_NAME ${GPU_GC_name}
      COMPILE_FLAGS "${GPU_GC_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(GCMC ws2_32)
   endif()
endif()

if(ENSEMBLE_GPU_NPT)
   add_executable(GPU_NPT ${sources} ${headers} ${libHeaders} ${cudaSources})
   set_target_properties(GPU_NPT PROPERTIES 
      OUTPUT_NAME ${GPU_NPT_name}
      COMPILE_FLAGS "${GPU_NPT_flags}")
   if(WIN32)
      #needed for hostname
      target_link_libraries(NPT ws2_32)
   endif()
endif()
