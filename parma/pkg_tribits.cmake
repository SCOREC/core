TRIBITS_PACKAGE(SCORECparma)

SET(PARMA_EXTERNAL_HEADERS parma.h)

SET(API_SOURCE parma.cc)

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
  diffMC/parma_elmBalancer.cc
  diffMC/parma_elmBdrySides.cc
  diffMC/parma_elmSideSides.cc
  diffMC/parma_edgeEqVtxSelector.cc
  diffMC/parma_ltSelector.cc
  diffMC/parma_elmLtVtxEdgeSelector.cc
  diffMC/parma_elmSelector.cc
  diffMC/parma_vtxSides.cc
  diffMC/parma_entWeights.cc
  diffMC/parma_ghost.cc
  diffMC/parma_ghostElement.cc
  diffMC/parma_ghostOwner.cc
  diffMC/parma_ghostWeights.cc
  diffMC/parma_ghostMPAS.cc
  diffMC/parma_ghostMPASWeights.cc
  diffMC/parma_vtxPtnWriter.cc
  diffMC/parma_graphDist.cc
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
  diffMC/parma_preserveTargets.cc
  diffMC/parma_vtxEdgeTargets.cc
  diffMC/parma_elmLtVtxEdgeTargets.cc
  diffMC/parma_vtxEdgeElmBalancer.cc
  diffMC/parma_vtxElmBalancer.cc
  diffMC/parma_elmLtVtxEdgeBalancer.cc
  diffMC/zeroOneKnapsack.c
  diffMC/maximalIndependentSet/misLuby.cc
  diffMC/maximalIndependentSet/mersenne_twister.cc
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

SET(parmaDepLibs ${APF_LIBS})

TRIBITS_ADD_LIBRARY(
  parma
  SOURCES ${DIFFMC_SOURCES} ${RIB_SOURCES} ${GROUP_SOURCES} ${API_SOURCE}
  HEADERS ${PARMA_EXTERNAL_HEADERS})

TRIBITS_PACKAGE_POSTPROCESS()
