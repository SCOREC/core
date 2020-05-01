tribits_package(SCORECem)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES
  emResidualFunctionals.cc)

set(HEADERS
  em.h)

#Library
tribits_add_library(
  em
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
