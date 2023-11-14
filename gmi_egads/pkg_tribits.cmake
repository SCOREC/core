tribits_package(SCORECgmi_egads)

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
