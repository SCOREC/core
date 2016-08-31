if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
set(Trilinos_PREFIX CACHE STRING "Trilinos installation directory")
option(ENABLE_STK "Enable the stk interface [ON|OFF]" OFF)
option(ENABLE_STK_MESH "Enable the stk mesh interface [ON|OFF]" OFF)
message(STATUS "ENABLE_STK: " ${ENABLE_STK})
message(STATUS "ENABLE_STK_MESH: " ${ENABLE_STK_MESH})

# Only install the package if enabled
if(NOT ENABLE_STK)
  return()
endif()

# Find the shards Trilinos package
set(HAS_STK OFF)
dep_macro(Shards "${Trilinos_PREFIX}")

# If enabling more mesh features, find additional Trilinos packages
if(ENABLE_STK_MESH)
  set(HAS_STK ON)
  dep_macro(STKIO "${Trilinos_PREFIX}")
  dep_macro(STKMesh "${Trilinos_PREFIX}")
  dep_macro(SEACASIoss "${Trilinos_PREFIX}")
endif()

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

if(HAS_STK)
  set(HEADERS ${HEADERS} apfSTK.h)
  set(SOURCES ${SOURCES} apfExodusOutput.cc)
endif()

# Add the apf_stk library
add_library(apf_stk ${SOURCES})

target_link_libraries(apf_stk
  PRIVATE gmi
  PUBLIC apf)
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
if(HAS_STK)
  target_link_libraries(apf_stk PUBLIC
    ${STKIO_LIBRARIES}
    ${STKMesh_LIBRARIES}
    ${STKMesh_TPL_LIBRARIES}
    ${SEACASIoss_LIBRARIES}
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

target_macro(apf_stk)

deps_end_macro()
