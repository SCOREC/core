# - Try to find parmetis
# Once done this will define
#  PARMETIS_FOUND - System has PARMETIS
#  PARMETIS_INCLUDE_DIRS - The PARMETIS include directories
#  PARMETIS_LIBRARIES - The libraries needed to use PARMETIS
#  PARMETIS_DEFINITIONS - Compiler switches required for using PARMETIS

set(PARMETIS_PREFIX "" CACHE STRING "ParMETIS install directory")

find_path(PARMETIS_INCLUDE_DIR parmetis.h PATHS "${PARMETIS_PREFIX}/include")
if(NOT EXISTS "${PARMETIS_INCLUDE_DIR}")
  message(FATAL_ERROR "parmetis include dir not found")
endif()

find_library(PARMETIS_LIBRARY parmetis PATHS "${PARMETIS_PREFIX}/lib")
if(NOT EXISTS "${PARMETIS_LIBRARY}")
  message(FATAL_ERROR "parmetis library not found")
endif()

find_library(METIS_LIBRARY metis PATHS "${PARMETIS_PREFIX}/lib")
if(NOT EXISTS "${METIS_LIBRARY}")
  message(FATAL_ERROR "metis library not found")
endif()

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PARMETIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PARMETIS  DEFAULT_MSG
    PARMETIS_LIBRARY METIS_LIBRARY PARMETIS_INCLUDE_DIR)

mark_as_advanced(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY METIS_LIBRARY)
