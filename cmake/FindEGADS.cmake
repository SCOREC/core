# - Try to find EGADS
# Once done this will define
#  EGADS_FOUND - System has EGADS
#  EGADS_INCLUDE_DIRS - The EGADS include directories
#  EGADS_LIBRARIES - The libraries needed to use EGADS
#  EGADS_DEFINITIONS - Compiler switches required for using EGADS

#set(EGADS_PREFIX "${EGADS_PREFIX_DEFAULT}" CACHE STRING "EGADS install directory")
#if(EGADS_PREFIX)
#  message(STATUS "EGADS_PREFIX: ${EGADS_PREFIX}")
#endif()

find_path(EGADS_INCLUDE_DIR egads.h PATHS "${EGADS_DIR}/include")

option(EGADS_LITE "Use EGADS_LITE" OFF)
if (EGADS_LITE)
  message(STATUS "Using EGADSlite")
  find_library(EGADS_LIBRARY NAMES egadslite PATHS "${EGADS_DIR}/lib")
else()
  message(STATUS "Using EGADS")
  find_library(EGADS_LIBRARY NAMES egads libegads.so libegads.dylib PATHS "${EGADS_DIR}/lib")
endif()

set(EGADS_LIBRARIES ${EGADS_LIBRARY} )
set(EGADS_INCLUDE_DIRS ${EGADS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EGADS_FOUND to TRUE
# if all listed variables are TRUE
message(STATUS "beep boop: ${EGADS_LIBRARY}")
find_package_handle_standard_args(
    EGADS
    DEFAULT_MSG
    EGADS_LIBRARY EGADS_INCLUDE_DIR
)

mark_as_advanced(EGADS_INCLUDE_DIR EGADS_LIBRARY)
