cmake_minimum_required(VERSION 2.8.6)

project(apf_mpas)
option(ENABLE_MPAS "Enable MPAS interface" OFF)

if (ENABLE_MPAS)
  find_package(parma PATHS ${CMAKE_BINARY_DIR})
  set(NETCDF_CXX "YES")
  find_package(NetCDF REQUIRED)
  find_package(HDF5 REQUIRED HL C)
  set(APF_MPAS_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${PARMA_INCLUDE_DIRS})
  set(DEP_LIBS ${PARMA_LIBS} 
      ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES})
  set(MPAS_LIBS apf_mpas ${DEP_LIBS})
endif()

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

if (NOT ENABLE_MPAS)
  return()
endif ()

include_directories(${APF_MPAS_INCLUDE_DIRS})
include_directories(${NETCDF_INCLUDES})
include_directories(${NETCDF_INCLUDES_CXX})

#Sources & Headers
set(SOURCES
   apfMPAS.cc)

#Library
add_library(apf_mpas ${SOURCES})
target_link_libraries(apf_mpas ${DEP_LIBS})

#Install
INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_mpas.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS  
   apf_mpas
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)
