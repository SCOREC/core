tribits_package(SCORECree)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES
  reeResidualFunctionals.cc
  reeFluxCorrection.cc
  reeCorrectedFlux.cc
  reeEstimateError.cc
  reeSizeField.cc)

set(HEADERS
  ree.h)

#Library
tribits_add_library(
  ree
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
