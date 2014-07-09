cmake_minimum_required(VERSION 2.8.6)

project(mpas)
find_package(Threads)

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/mpasConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/mpasConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "0.1")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libmpas.pc.in"
    "${PROJECT_BINARY_DIR}/libmpas.pc"
    @ONLY)
#pkgconfig }

#directory containing pcu header files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

#Sources & Headers
set(SOURCES
   mpas.cc)

set(HEADERS
   mpas.h)

#Library
add_library(mpas ${SOURCES})

#Install
INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libmpas.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS  
   mpas
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)

