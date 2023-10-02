tribits_package(SCORECapf_stk)

# determine if configured with stk
SET(ENABLE_STK_MESH OFF)
if(SCORECapf_stk_ENABLE_STKIO AND
   SCORECapf_stk_ENABLE_STKMesh AND
   SCORECapf_stk_ENABLE_SEACASIoss)
  SET(ENABLE_STK_MESH ON)
endif()

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

if(ENABLE_STK_MESH)
  set(HEADERS ${HEADERS} apfSTK.h)
  set(SOURCES ${SOURCES} apfExodusOutput.cc)
endif()

#Library
tribits_add_library(
  apf_stk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
