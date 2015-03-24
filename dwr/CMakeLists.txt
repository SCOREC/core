cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros) 
  include(TribitsLibraryMacros) 
  tribits_subpackage(dwr) 

  SET(SOURCES )
  SET(HEADERS )

  # dwr dirctory
  SET(DWR_DIR ${CMAKE_CURRENT_SOURCE_DIR})

  # include directories
  SET(INCLUDE_DIRS
    ${DWR_DIR}
    ${DWR_DIR}/qoi
    ${DWR_DIR}/elasticity)
  include_directories(${INCLUDE_DIRS})

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/dwrConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/dwrConfig.cmake")

  #pkgconfig {
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
  set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
  set(includedir "${CMAKE_INSTALL_PREFIX}/include")
  set(PACKAGE_VERSION "1.0")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libdwr.pc.in"
    "${PROJECT_BINARY_DIR}/libdwr.pc"
    @ONLY)
  #pkgconfig }

  SET(BASE_SOURCES
    dwr.cc
    dwrUtils.cc
    dwrHierarchic.cc
    dwrLinearSystem.cc)

  SET(QOI_SOURCES
    qoi/dwrVectorL2QOI.cc)

  SET(ELASTICITY_SOURCES
    elasticity/dwrElasticityLHS.cc
    elasticity/dwrElasticityErr.cc
    elasticity/dwrElasticityProblem.cc)
    
  set(SOURCES
    ${BASE_SOURCES}
    ${QOI_SOURCES}
    ${ELASTICITY_SOURCES})
    
  set(HEADERS
    dwr.h
    dwrHierarchic.h
    elasticity/dwrElasticityProblem.h)

  #Library
  tribits_add_library(
    dwr
    HEADERS ${HEADERS}
    SOURCES ${SOURCES})

  tribits_subpackage_postprocess()
endif(BUILD_IN_TRILINOS) 
