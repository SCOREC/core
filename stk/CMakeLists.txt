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
find_package(Shards REQUIRED PATHS ${Trilinos_PREFIX})

# If enabling more mesh features, find additional Trilinos packages
if(ENABLE_STK_MESH)
  set(HAS_STK ON)
  find_package(STKIO REQUIRED PATHS ${Trilinos_PREFIX})
  find_package(STKMesh REQUIRED PATHS ${Trilinos_PREFIX})
  find_package(STKIoss REQUIRED PATHS ${Trilinos_PREFIX})
endif()

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
set(
  apfAlbany.h
)

# Set the link libraries and include directories
set(APF_STK_INCS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${Shards_INCLUDE_DIRS}
)
set(APF_STK_LIBS
  apf
  ${Shards_LIBRARIES}
)

# Do some extra work if has stk is on
if(HAS_STK)
  set(HEADERS ${HEADERS} apfSTK.h)
  set(SOURCES ${SOURCES} apfExodusOutput.cc)
  set(APF_STK_LIBS
    ${APF_STK_LIBS}
    ${STKIO_LIBRARIES}
    ${STKMesh_LIBRARIES}
    ${STKIoss_LIBRARIES}
  )
  set(APF_STK_INCS
    ${APF_STK_INCS}
    ${STKIO_INCLUDE_DIRS}
    ${STKMesh_INCLUDE_DIRS}
    ${STKIoss_INCLUDE_DIRS}
  )
endif()

# Add the apf_stk library
add_library(apf_stk ${SOURCES})

# Inlcude directories
include_directories(${APF_STK_INCS})

# Link this package to these libraries
target_link_libraries(apf_stk ${APF_STK_LIBS})

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET apf_stk
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS apf_stk EXPORT apf_stk-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT apf_stk-targets DESTINATION lib/scorec)
