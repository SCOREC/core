cmake_minimum_required(VERSION 2.8.6)

IF(DEFINED TRIBITS_PACKAGE)
 SET(BUILD_IN_TRILINOS ON)
ENDIF()

if(BUILD_IN_TRILINOS)
   # determine if Trilinos has been configured with stk
   SET(HAS_STK OFF)
   if(Trilinos_ENABLE_STKIO AND
      Trilinos_ENABLE_STKMesh AND
      Trilinos_ENABLE_SEACASIoss)
     SET(HAS_STK ON)
   endif()

   tribits_package(SCORECapf_stk) 
else()
   project(apf_stk)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  set(APF_STK_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${APF_INCLUDE_DIRS}
      ${STK_INCLUDE_DIRS})
  set(APF_STK_LIBS
      apf_stk
      ${APF_LIBS}
      ${STK_LIBRARIES})
  include_directories(${APF_STK_INCLUDE_DIRS})
else(NOT BUILD_IN_TRILINOS)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_stkConfig.h.in"
    "${PROJECT_BINARY_DIR}/apf_stkConfig.h")
include_directories("${PROJECT_BINARY_DIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_stkConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apf_stkConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libapf_stk.pc.in"
    "${PROJECT_BINARY_DIR}/libapf_stk.pc"
    @ONLY)
#pkgconfig }

#Sources & Headers
set(SOURCES
  apfMeshSTK.cc
  apfSTK.cc)

set(HEADERS apfAlbany.h)
if(HAS_STK)
  set(HEADERS ${HEADERS} apfSTK.h)
  set(SOURCES ${SOURCES} apfExodusOutput.cc)
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
     apf_stk
     HEADERS ${HEADERS}
     SOURCES ${SOURCES})
else()
   add_library(apf_stk ${SOURCES})
   target_link_libraries(apf_stk ${APF_STK_LIBS})
endif()

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_stk.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      apf_stk
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_package_postprocess()
endif()

