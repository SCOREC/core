if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

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
  maRegionCollapse.cc
  maMatchedCollapse.cc
  maLayerCollapse.cc
  maMatch.cc
  maSolutionTransfer.cc
  maSolutionTransferHelper.cc
  maSnap.cc
  maEdgeSwap.cc
  maShape.cc
  maShapeHandler.cc
  maQuality.cc
  maSplits.cc
  maFaceSplit.cc
  maDoubleSplitCollapse.cc
  maSingleSplitCollapse.cc
  maFaceSplitCollapse.cc
  maShortEdgeRemover.cc
  maVertRemover.cc
  maSnapper.cc
  maMatchedSnapper.cc
  maBalance.cc
  maLayer.cc
  maCrawler.cc
  maTetrahedronize.cc
  maLayerSnap.cc
  maMap.cc
  maReposition.cc
  maExtrude.cc
  maDBG.cc
  maStats.cc
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
  maDBG.h
  maStats.h
)

# Add the ma library
add_library(ma ${SOURCES})

# Include directories
target_include_directories(ma INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

# Link this library to these libraries
target_link_libraries(ma
   PUBLIC
     apf
     mds
     parma
     apf_cap
     apf_zoltan
     apf_metis
     pcu
   )

scorec_export_library(ma)

bob_end_subdir()
