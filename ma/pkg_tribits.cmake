tribits_package(SCORECma)

#Sources & Headers
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
  maFixShape.cc
  maFaceSwap.cc
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

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Library
tribits_add_library(
   ma
   HEADERS ${HEADERS}
   SOURCES ${SOURCES})

TRIBITS_PACKAGE_POSTPROCESS()
