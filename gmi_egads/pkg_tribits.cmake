tribits_package(SCORECgmi_egads)

# these ENABLE vars are made by tribits based on ./cmake/Dependencies.cmake
option(EGADS_LITE "Enable EGADSlite" OFF)
# this file brings EGADS_LITE from CMake to C++
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gmi_egads_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/gmi_egads_config.h")
# make sure the compiler can find the above header
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES gmi_egads.c)
set(HEADERS gmi_egads.h)

#Library
tribits_add_library(
  gmi_egads
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
