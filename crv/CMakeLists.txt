cmake_minimum_required(VERSION 2.8.6)

project(crv)

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(ma PATHS ${CMAKE_BINARY_DIR})
  
  set(CRV_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${MA_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS} ${MA_LIBS})
  set(CRV_LIBS crv ${DEP_LIBS})
endif()

#Sources & Headers
set(CRV_SOURCES
  crv.cc
  crvBernstein.cc
  crvBezier.cc
  crvBezierPoints.cc
  crvBlended.cc
  crvCurveMesh.cc
  crvSnap.cc
  crvSubdivision.cc
  crvTables.cc
  crvQuality.cc
  crvVtk.cc)

set(CRV_HEADERS
  crv.h)

include_directories(${CRV_INCLUDE_DIRS})

#Library
add_library(crv ${CRV_SOURCES})
target_link_libraries(crv ${DEP_LIBS})

INSTALL(FILES ${CRV_HEADERS} DESTINATION include)
INSTALL(TARGETS 
  crv
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/crvConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/crvConfig.cmake")

