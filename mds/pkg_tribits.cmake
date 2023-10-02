tribits_package(SCORECmds)

set(MDS_SET_MAX 256 CACHE STRING "Buffer size for adjacency computation")
set(MDS_ID_TYPE "int" CACHE STRING "Interal identifier integer type")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mds_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/mds_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#Sources & Headers
set(MDS_SOURCES
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
  mdsUgrid.cc)

set(MDS_HEADERS
  apfMDS.h
  apfBox.h
)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Library
tribits_add_library(
   mds
   HEADERS ${MDS_HEADERS}
   SOURCES ${MDS_SOURCES})

tribits_package_postprocess()
