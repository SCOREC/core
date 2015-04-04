cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(dsp) 
else()
   project(dsp)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(ma PATHS ${CMAKE_BINARY_DIR})
  set(DSP_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MA_INCLUDE_DIRS})
  set(DEP_LIBS ${MA_LIBS})
  set(DSP_LIBS dsp ${DEP_LIBS})
endif()

#Sources & Headers
set(SOURCES
  dspAdapters.cc
  dspSmoothers.cc
  dsp.cc)

set(HEADERS
  dspAdapters.h
  dspSmoothers.h
  dsp.h)

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else()
  include_directories(${DSP_INCLUDE_DIRS})
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      dsp
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(dsp ${SOURCES})
   target_link_libraries(dsp ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libdsp.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS 
     dsp
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/dspConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/dspConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libdsp.pc.in"
    "${PROJECT_BINARY_DIR}/libdsp.pc"
    @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
   TRIBITS_SUBPACKAGE_POSTPROCESS()
endif()
