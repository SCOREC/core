if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

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

if(ENABLE_CGNS)
  set(SOURCES ${SOURCES} mdsCGNS.cc)
endif(ENABLE_CGNS)

# Package headers
set(HEADERS
  apfMDS.h
  apfBox.h
)

# Add the mds library
add_library(mds ${SOURCES})

# Include directories
target_include_directories(mds INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )
# for the generated mds_config.h
# private because no one outside MDS needs it
target_include_directories(mds PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    )

# Link this library to these libraries
target_link_libraries(mds
   PUBLIC
     pcu
     gmi
     apf
   )

if(ENABLE_CGNS)
  message(STATUS ${CGNS_LIBRARIES})
  target_link_libraries(mds PRIVATE ${CGNS_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif(ENABLE_CGNS)

scorec_export_library(mds)

bob_end_subdir()
