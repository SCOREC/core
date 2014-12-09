cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(gmi) 
else()
   project(gmi)
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/gmiConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/gmiConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "2.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libgmi.pc.in"
    "${PROJECT_BINARY_DIR}/libgmi.pc"
    @ONLY)
#pkgconfig }

#Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

#directory containing gmi header files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCES
   gmi.c
   agm.c
   gmi_base.c
   gmi_file.c
   gmi_lookup.c
   gmi_mesh.c
   gmi_null.c
   gmi_analytic.c)

set(HEADERS
   gmi.h
   agm.h
   gmi_base.h
   gmi_mesh.h
   gmi_null.h
   gmi_analytic.h)

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      gmi 
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(gmi ${SOURCES})
endif()

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libgmi.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      gmi
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif()
