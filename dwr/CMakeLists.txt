cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros) 
  include(TribitsLibraryMacros) 
  tribits_subpackage(dwr) 

  # THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

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

  #Base Sources
  set(SOURCES
    dwr.cc
    dwrHierarchic.cc
    dwrLinearSystem.cc)

  #Header
  set(HEADERS
    dwr.h
    dwrHierarchic.h)

  #Library
  tribits_add_library(
    dwr
    HEADERS ${HEADERS}
    SOURCES ${SOURCES})

  tribits_subpackage_postprocess()
endif(BUILD_IN_TRILINOS) 
