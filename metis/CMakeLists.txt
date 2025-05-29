if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package sources
set(SOURCES
  apfMETIS.cc
  apfMETIScommon.cc
  apfMETISsplitter.cc
  apfMETISbalancer.cc
)

# Package headers
set(HEADERS
  apfMETIS.h
  apfMETIScommon.h
  apfMETISsplitter.h
  apfMETISbalancer.h
)

# Add the apf_metis library
add_library(apf_metis ${SOURCES})


# Include directories
target_include_directories(apf_metis INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Link this library to these libraries
target_link_libraries(apf_metis PUBLIC apf)
target_link_libraries(apf_metis PRIVATE apf_zoltan) # FIXME: tagOpposites

set(METIS_PREFIX "" CACHE PATH "Prefix for METIS")
if(METIS_PREFIX)
  target_include_directories(apf_metis PUBLIC "${METIS_PREFIX}/include")
  target_link_directories(apf_metis INTERFACE "${METIS_PREFIX}/lib")
  target_link_libraries(apf_metis INTERFACE metis)
endif()

scorec_export_library(apf_metis)

bob_end_subdir()
