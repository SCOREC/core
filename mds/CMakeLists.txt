# Package options
set(MDS_SET_MAX "256" CACHE STRING "Buffer size for adjacency computation")
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

# Include directories
target_include_directories(apf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
# for the generated mds_config.h
# private because no one outside MDS needs it
target_include_directories(apf PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    )

# Link this library to these libraries
target_link_libraries(mds
# these headers are visible from public MDS headers
   PRIVATE pcu
# these headers are visible from public MDS headers
   PUBLIC gmi apf
   )

# Tell CMake to add this package to the SCOREC project
install(TARGETS mds EXPORT mds-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT mds-targets DESTINATION lib/scorec)
