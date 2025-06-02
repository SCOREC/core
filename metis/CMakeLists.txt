if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

find_package(METIS)

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
target_compile_features(apf_metis PRIVATE cxx_std_14)
target_include_directories(apf_metis INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(apf_metis PUBLIC apf)
target_link_libraries(apf_metis PRIVATE apf_zoltan) # FIXME: tagOpposites
if(METIS_FOUND)
  target_compile_definitions(apf_metis PUBLIC PUMI_HAS_METIS)
  target_link_libraries(apf_metis PRIVATE METIS::METIS)
endif()

scorec_export_library(apf_metis)
bob_end_subdir()
