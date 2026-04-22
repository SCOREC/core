tribits_package(SCORECapf_sim)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES apfSIM.cc)
set(HEADERS apfSIM.h)

#Library
tribits_add_library(
  apf_sim
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
