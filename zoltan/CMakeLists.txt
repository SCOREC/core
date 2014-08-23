cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros) 
  include(TribitsLibraryMacros) 
  tribits_subpackage(apf_zoltan) 
  set(ENABLE_ZOLTAN ON)
else()
  project(apf_zoltan)
  option(ENABLE_ZOLTAN "Enable Zoltan interface" OFF)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  if(ENABLE_ZOLTAN)
    find_package(Zoltan QUIET REQUIRED)
  endif()
  set(APF_ZOLTAN_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${APF_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS})
  if(ENABLE_ZOLTAN)
    set(DEP_LIBS ${DEP_LIBS} ${ZOLTAN_LIBRARIES} ${PARMETIS_LIBRARIES})
  endif()
  set(APF_ZOLTAN_LIBS apf_zoltan ${DEP_LIBS})
  include_directories(${APF_ZOLTAN_INCLUDE_DIRS})
  if(ENABLE_ZOLTAN)
    include_directories(${ZOLTAN_INCLUDE_DIRS})
    include_directories(${PARMETIS_INCLUDE_DIRS})
  endif()
else(NOT BUILD_IN_TRILINOS)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
endif()


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_zoltanConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apf_zoltanConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libapf_zoltan.pc.in"
    "${PROJECT_BINARY_DIR}/libapf_zoltan.pc"
    @ONLY)
#pkgconfig }

#Sources & Headers
if(ENABLE_ZOLTAN)
  set(SOURCES
    apfInterElement.cc
    apfZoltan.cc
    apfZoltanMesh.cc
    apfZoltanCallbacks.cc)
else()
  set(SOURCES
    apfInterElement.cc
    apfZoltanEmpty.cc)
endif()

set(HEADERS
  apfZoltan.h)

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
     apf_zoltan
     HEADERS ${HEADERS}
     SOURCES ${SOURCES})
else()
   add_library(apf_zoltan ${SOURCES})
   target_link_libraries(apf_zoltan ${DEP_LIBS})
endif()

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_zoltan.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      apf_zoltan
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif()

