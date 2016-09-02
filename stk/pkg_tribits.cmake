# determine if Trilinos has been configured with stk
SET(ENABLE_STK_MESH OFF)
if(Trilinos_ENABLE_STKIO AND
   Trilinos_ENABLE_STKMesh AND
   Trilinos_ENABLE_SEACASIoss)
  SET(ENABLE_STK_MESH ON)
endif()

tribits_package(SCORECapf_stk) 

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_stkConfig.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/apf_stkConfig.h")
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

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
tribits_add_library(
  apf_stk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
