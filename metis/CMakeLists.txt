if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

option(PUMI_ENABLE_METIS "Enable METIS interface" ${ENABLE_METIS})
message(VERBOSE "PUMI_ENABLE_METIS: ${PUMI_ENABLE_METIS}")
xsdk_add_tpl(METIS)
if(PUMI_ENABLE_METIS)
  find_package(METIS REQUIRED)
endif()

# Package sources
if(PUMI_ENABLE_METIS)
  set(SOURCES
    apfMETIS.cc
    apfMETIScommon.cc
    apfMETISsplitter.cc
    apfMETISbalancer.cc
  )
else()
  set(SOURCES apfMETISempty.cc)
endif()

# Package headers
set(HEADERS apfMETIS.h)

# Add the apf_metis library
add_library(apf_metis ${SOURCES})
target_include_directories(apf_metis INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(apf_metis PUBLIC apf)
if(PUMI_ENABLE_METIS)
  target_compile_definitions(apf_metis PUBLIC PUMI_HAS_METIS)
  target_link_libraries(apf_metis PRIVATE "${METIS_LIBRARIES}")
  target_include_directories(apf_metis PRIVATE "${METIS_INCLUDE_DIRS}")
endif()

scorec_export_library(apf_metis)
bob_end_subdir()
