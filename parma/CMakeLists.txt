cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  INCLUDE(TribitsSubPackageMacros) 
  INCLUDE(TribitsLibraryMacros) 
  TRIBITS_SUBPACKAGE(parma) 
else()
  project(parma)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
endif()

SET(PARMA_EXTERNAL_HEADERS parma.h)

SET(API_SOURCE parma.cc)

######################################################################
SET(DIFFMC_SOURCES  
  diffMC/parma_balancer.cc
  diffMC/parma_bdryVtx.cc
  diffMC/parma_centroidDiffuser.cc
  diffMC/parma_centroids.cc
  diffMC/parma_centroidSelector.cc
  diffMC/parma_commons.cc
  diffMC/parma_components.cc
  diffMC/parma_dcpart.cc
  diffMC/parma_dcpartFixer.cc
  diffMC/parma_dijkstra.cc
  diffMC/parma_edgeBalancer.cc
  diffMC/parma_edgeSelector.cc
  diffMC/parma_elmBalancer.cc
  diffMC/parma_elmBdrySides.cc
  diffMC/parma_elmSideSides.cc
  diffMC/parma_elmLtVtxSelector.cc
  diffMC/parma_elmSelector.cc
  diffMC/parma_vtxSides.cc
  diffMC/parma_entWeights.cc
  diffMC/parma_ghost.cc
  diffMC/parma_ghostOwner.cc
  diffMC/parma_ghostWeights.cc
  diffMC/parma_graphDist.cc
  diffMC/parma_hpsBalancer.cc
  diffMC/parma_monitor.cc
  diffMC/parma_sides.cc
  diffMC/parma_step.cc
  diffMC/parma_stop.cc
  diffMC/parma_shapeOptimizer.cc
  diffMC/parma_shapeTargets.cc
  diffMC/parma_shapeSelector.cc
  diffMC/parma_vtxBalancer.cc
  diffMC/parma_vtxSelector.cc
  diffMC/parma_weightTargets.cc
  diffMC/parma_weightSideTargets.cc
  diffMC/parma_vtxElmTargets.cc
  diffMC/parma_vtxEdgeElmTargets.cc
  diffMC/parma_vtxEdgeElmBalancer.cc
  diffMC/parma_vtxElmBalancer.cc
  diffMC/zeroOneKnapsack.c
  diffMC/maximalIndependentSet/misLuby.cc
  )

SET(RIB_SOURCES  
  rib/parma_rib.cc 
  rib/parma_mesh_rib.cc 
  )

SET(GROUP_SOURCES  
  group/parma_group.cc 
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(diffMC)

if(NOT BUILD_IN_TRILINOS) 
  set(parmaIncDirs 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    )
  INCLUDE_DIRECTORIES(${parmaIncDirs}
    ${PUMI_GEOM_MESHMODEL_INCLUDE_DIRS})
endif()

SET(parmaDepLibs ${APF_LIBS})

if(BUILD_IN_TRILINOS)
  TRIBITS_ADD_LIBRARY(
    parma
    SOURCES ${DIFFMC_SOURCES} ${RIB_SOURCES} ${GROUP_SOURCES} ${API_SOURCE}
    HEADERS ${PARMA_EXTERNAL_HEADERS})
else(BUILD_IN_TRILINOS)
  ADD_LIBRARY(parma 
    ${DIFFMC_SOURCES}
    ${RIB_SOURCES}
    ${GROUP_SOURCES}
    ${API_SOURCE})
  target_link_libraries(parma ${parmaDepLibs})
endif(BUILD_IN_TRILINOS)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/parmaConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/parmaConfig.cmake"
)

if(NOT BUILD_IN_TRILINOS)
  #pkgconfig {
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
  set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
  set(includedir "${CMAKE_INSTALL_PREFIX}/include")
  set(PACKAGE_VERSION "1.0")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libparma.pc.in"
    "${PROJECT_BINARY_DIR}/libparma.pc"
    @ONLY
    )
  #pkgconfig }
endif()

if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES ${PARMA_EXTERNAL_HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libparma.pc" DESTINATION lib/pkgconfig)
  INSTALL(TARGETS 
    parma
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
  TRIBITS_SUBPACKAGE_POSTPROCESS()
endif()
