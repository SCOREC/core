# Package sources
set(SOURCES
  ma.cc
  maInput.cc
  maAdapt.cc
  maMesh.cc
  maRefine.cc
  maLayerRefine.cc
  maLayerCoarsen.cc
  maTables.cc
  maLayerTables.cc
  maTemplates.cc
  maLayerTemplates.cc
  maCoarsen.cc
  maSize.cc
  maOperator.cc
  maCollapse.cc
  maMatchedCollapse.cc
  maLayerCollapse.cc
  maMatch.cc
  maSolutionTransfer.cc
  maSnap.cc
  maEdgeSwap.cc
  maShape.cc
  maShapeHandler.cc
  maQuality.cc
  maSplits.cc
  maDoubleSplitCollapse.cc
  maShortEdgeRemover.cc
  maVertRemover.cc
  maSnapper.cc
  maBalance.cc
  maLayer.cc
  maCrawler.cc
  maTetrahedronize.cc
  maLayerSnap.cc
  maMap.cc
  maReposition.cc
  maExtrude.cc
)

# Package headers
set(HEADERS
  ma.h
  maInput.h
  maMesh.h
  maSize.h
  maShape.h
  maTables.h
  maSolutionTransfer.h
  maExtrude.h
)

# Add the ma library
add_library(ma ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link this package to these libraries
target_link_libraries(ma apf mds parma apf_zoltan)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET ma
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS ma EXPORT ma-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT ma-targets DESTINATION lib/scorec)
