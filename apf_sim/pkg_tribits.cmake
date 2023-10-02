tribits_package(SCORECapf_sim)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(SCORECapf_sim_ENABLE_SimField)
  add_definitions(-DUSE_FIELDSIM)
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.h.in"
    "${PROJECT_BINARY_DIR}/apf_simConfig.h")
#directory containing apf_simConfig.h
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

#Sources & Headers
set(SOURCES apfSIM.cc)
set(HEADERS apfSIM.h)

#Library
tribits_add_library(
  apf_sim 
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
