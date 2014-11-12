cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(spr) 
else()
   project(spr)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  set(SPR_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${APF_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS})
  set(SPR_LIBS spr ${DEP_LIBS})
endif()

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else(BUILD_IN_TRILINOS)
  include_directories(${SPR_INCLUDE_DIRS})
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/sprConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/sprConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libspr.pc.in"
    "${PROJECT_BINARY_DIR}/libspr.pc"
    @ONLY)
#pkgconfig }

#Sources & Headers
set(SOURCES
  sprGetGradIPField.cc
  sprSVD.cc
  sprQR.cc
  sprRecoverField.cc
  sprEstimateError.cc)

set(HEADERS
  spr.h)

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
     spr
     HEADERS ${HEADERS}
     SOURCES ${SOURCES})
else()
   add_library(spr ${SOURCES})
   target_link_libraries(spr ${DEP_LIBS})
endif()

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libspr.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      spr
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif()

