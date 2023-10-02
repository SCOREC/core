tribits_package(SCORECpumi)

#Sources & Headers
set(SOURCES
  GenTag.cc
  mPartEntityContainer.cc
  pumi_geom.cc
  pumi_gentity.cc
  pumi_ghost.cc
  pumi_gtag.cc
  pumi_mesh.cc
  pumi_mentity.cc
  pumi_mtag.cc
  pumi_sys.cc
)

set(HEADERS
  pumi_errorcode.h
  GenTag.h
  GenIterator.h
  pumi_list.h
  mPartEntityContainer.h
  pumi.h
)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Library
tribits_add_library(
   pumi
   HEADERS ${HEADERS}
   SOURCES ${SOURCES})

tribits_package_postprocess()
