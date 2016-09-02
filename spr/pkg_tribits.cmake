tribits_package(SCORECspr)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES
  sprGetGradIPField.cc
  sprRecoverField.cc
  sprEstimateError.cc
  sprEstimateTargetError.cc)

set(HEADERS
  spr.h)

#Library
tribits_add_library(
  spr
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
