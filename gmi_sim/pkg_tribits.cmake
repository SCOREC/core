tribits_package(SCORECgmi_sim)

# these ENABLE vars are made by tribits based on ./cmake/Dependencies.cmake
set(SIM_PARASOLID ${SCORECgmi_sim_ENABLE_SimParasolid})
set(SIM_ACIS ${SCORECgmi_sim_ENABLE_SimAcis})
# this file brings SIM_PARASOLID and SIM_ACIS from CMake to C++
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gmi_sim_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/gmi_sim_config.h")
# make sure the compiler can find the above header
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES gmi_sim.cc)
set(HEADERS gmi_sim.h)

#Library
tribits_add_library(
  gmi_sim 
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
