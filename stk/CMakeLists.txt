if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(ENABLE_STK "Enable the stk interface [ON|OFF]" OFF)
option(ENABLE_STK_MESH "Enable the stk mesh interface [ON|OFF]" OFF)
message(STATUS "ENABLE_STK: " ${ENABLE_STK})
message(STATUS "ENABLE_STK_MESH: " ${ENABLE_STK_MESH})

# Only install the package if enabled
if(NOT ENABLE_STK)
  return()
endif()

# Find the shards Trilinos package
set(SCOREC_USE_Shards_DEFAULT ${ENABLE_STK})
set(SCOREC_USE_STKIO_DEFAULT ${ENABLE_STK_MESH})
set(SCOREC_USE_STKMesh_DEFAULT ${ENABLE_STK_MESH})
set(SCOREC_USE_SEACASIoss_DEFAULT ${ENABLE_STK_MESH})
set(Shards_PREFIX_DEFAULT "${Trilinos_PREFIX}")
set(STKIO_PREFIX_DEFAULT "${Trilinos_PREFIX}")
set(STKMesh_PREFIX_DEFAULT "${Trilinos_PREFIX}")
set(SEACASIoss_PREFIX_DEFAULT "${Trilinos_PREFIX}")
bob_public_dep(Shards)
bob_public_dep(STKIO)
bob_public_dep(STKMesh)
bob_public_dep(SEACASIoss)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Configure the input header
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/apf_stkConfig.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/apf_stkConfig.h
)

# Package sources
set(SOURCES
  apfMeshSTK.cc
  apfSTK.cc
)

# Package headers
set(HEADERS
  apfAlbany.h
)

if(ENABLE_STK_MESH)
  set(HEADERS ${HEADERS} apfSTK.h)
  set(SOURCES ${SOURCES} apfExodusOutput.cc)
endif()

# Add the apf_stk library
add_library(apf_stk ${SOURCES})

target_link_libraries(apf_stk PUBLIC apf gmi)
target_link_libraries(apf_stk PUBLIC ${Shards_LIBRARIES})

# Include directories
target_include_directories(apf_stk INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
target_include_directories(apf_stk PUBLIC
    ${Shards_INCLUDE_DIRS}
    )
# Find the config header while building
target_include_directories(apf_stk PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    )

# Do some extra work if has stk is on
if(ENABLE_STK_MESH)
  target_link_libraries(apf_stk PUBLIC
    stk_io
    stk_io_util
    stk_mesh_base
    ${STKMesh_TPL_LIBRARIES}
    Ioss
    Ioex
  )
  target_include_directories(apf_stk
    SYSTEM
    PUBLIC
    ${STKIO_INCLUDE_DIRS}
    ${STKMesh_INCLUDE_DIRS}
    ${STKMesh_TPL_INCLUDE_DIRS}
    ${SEACASIoss_INCLUDE_DIRS}
  )
endif()

scorec_export_library(apf_stk)

bob_end_subdir()
