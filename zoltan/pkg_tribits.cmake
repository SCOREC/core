tribits_package(SCORECapf_zoltan)
set(ENABLE_ZOLTAN ON)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

#Sources & Headers
if(ENABLE_ZOLTAN)
  set(SOURCES
    apfInterElement.cc
    apfZoltan.cc
    apfZoltanMesh.cc
    apfZoltanCallbacks.cc)
else()
  set(SOURCES
    apfInterElement.cc
    apfZoltanEmpty.cc)
endif()

set(HEADERS
  apfZoltan.h)

#Library
tribits_add_library(
  apf_zoltan
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
