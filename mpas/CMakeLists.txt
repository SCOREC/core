cmake_minimum_required(VERSION 2.8.6)

project(apf_mpas)
option(ENABLE_MPAS "Enable MPAS interface" OFF)

if (NOT ENABLE_MPAS)
  return()
endif()

find_package(apf PATHS ${CMAKE_BINARY_DIR})
set(NETCDF_CXX "YES")
find_package(NetCDF REQUIRED)

set(MPAS_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS})
set(DEP_LIBS ${APF_LIBS} ${NETCDF_LIBRARIES})
set(MPAS_LIBS apf_mpas ${DEP_LIBS})

include_directories(${MPAS_INCLUDE_DIRS})
include_directories(${NETCDF_INCLUDES})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_mpasConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apf_mpasConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "0.1")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libapf_mpas.pc.in"
    "${PROJECT_BINARY_DIR}/libapf_mpas.pc"
    @ONLY)
#pkgconfig }

#Sources & Headers
set(SOURCES
   apfMPAS.cc)

set(HEADERS
   apfMPAS.h)

#Library
add_library(apf_mpas ${SOURCES})

#Install
INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_mpas.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS  
   apf_mpas
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)

