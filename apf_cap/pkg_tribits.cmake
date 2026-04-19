tribits_package(SCORECapf_cap)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#directory containing simModSuiteConfig.h
include_directories("${PROJECT_BINARY_DIR}")

#Sources & Headers
set(SOURCES apfCAP.cc)
set(HEADERS apfCAP.h)

#Library
tribits_add_library(
  apf_cap 
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
