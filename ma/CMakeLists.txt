cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(ma) 
else()
   project(ma)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(parma PATHS ${CMAKE_BINARY_DIR})
  find_package(apf_zoltan PATHS ${CMAKE_BINARY_DIR})
  set(MA_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${PARMA_INCLUDE_DIRS}
    ${APF_ZOLTAN_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS} ${PARMA_LIBS} ${APF_ZOLTAN_LIBS})
  set(MA_LIBS ma ${DEP_LIBS})
endif()

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
  maMap.cc)

set(HEADERS
  ma.h
  maInput.h
  maMesh.h
  maSize.h
  maSolutionTransfer.h)

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else()
  include_directories(${MA_INCLUDE_DIRS})
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      ma
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(ma ${SOURCES})
   target_link_libraries(ma ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libma.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS 
     ma
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/maConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/maConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libma.pc.in"
    "${PROJECT_BINARY_DIR}/libma.pc"
    @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
   TRIBITS_SUBPACKAGE_POSTPROCESS()
endif()
