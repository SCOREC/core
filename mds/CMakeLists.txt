# Package options
set(MDS_SET_MAX 256 CACHE STRING "Buffer size for adjacency computation")
set(MDS_ID_TYPE "int" CACHE STRING "Interal identifier integer type")
message(STATUS "MDS_SET_MAX: ${MDS_SET_MAX}")
message(STATUS "MDS_ID_TYPE: ${MDS_ID_TYPE}")

# Configure the input header
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/mds_config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/mds_config.h
)

# Package sources
set(SOURCES
  mds.c
  mds_apf.c
  mds_net.c
  mds_order.c
  mds_smb.c
  mds_tag.c
  apfMDS.cc
  apfPM.cc
  apfBox.cc
  mdsANSYS.cc
  mdsGmsh.cc
  mdsUgrid.cc
)

# Package headers
set(HEADERS
  apfMDS.h
  apfBox.h
)

# Add the mds library
add_library(mds ${SOURCES})

# Inlcude directories
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu
  ${CMAKE_CURRENT_SOURCE_DIR}/../gmi
  ${CMAKE_CURRENT_SOURCE_DIR}/../apf
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu/reel
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu/noto
)

# Link this package to these libraries
target_link_libraries(mds pcu gmi apf)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET mds
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS mds EXPORT mds-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT mds-targets DESTINATION lib/scorec)
