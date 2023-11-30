tribits_package(SCORECgmi_cap)

# make sure the compiler can find the above header
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Sources & Headers
set(SOURCES gmi_cap.cc)
set(HEADERS gmi_cap.h)

#Library
tribits_add_library(
  gmi_cap 
  HEADERS ${HEADERS}
  SOURCES ${SOURCES})

tribits_package_postprocess()
