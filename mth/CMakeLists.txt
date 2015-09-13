cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

if(DEFINED TRIBITS_PACKAGE)
 set(BUILD_IN_TRILINOS ON)
endif()

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros)
  include(TribitsLibraryMacros)
  tribits_subpackage(Can) 
else()
   project(mth)
endif()

set(MTH_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MTH_LIBS)

set(MTH_HEADERS
  mth.h
  mthVector.h
  mthMatrix.h
)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${MTH_INCLUDE_DIRS})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/mthConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/mthConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libmth.pc.in"
    "${PROJECT_BINARY_DIR}/libmth.pc"
    @ONLY)
#pkgconfig }

if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES ${MTH_HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libmth.pc" DESTINATION lib/pkgconfig)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif(BUILD_IN_TRILINOS)
