# - Try to find zoltan
# Once done this will define
#  ZOLTAN_FOUND - System has ZOLTAN
#  ZOLTAN_INCLUDE_DIRS - The ZOLTAN include directories
#  ZOLTAN_LIBRARIES - The libraries needed to use ZOLTAN
#  ZOLTAN_DEFINITIONS - Compiler switches required for using ZOLTAN

set(ZOLTAN_PREFIX "${ZOLTAN_PREFIX_DEFAULT}" CACHE STRING "Zoltan install directory")
if(ZOLTAN_PREFIX)
  message(STATUS "ZOLTAN_PREFIX ${ZOLTAN_PREFIX}")
endif()

if(DEFINED ZOLTAN_PREFIX)
  find_path(ZOLTAN_INCLUDE_DIR zoltan.h
    PATHS "${ZOLTAN_PREFIX}/include"
    NO_DEFAULT_PATH
  )
  find_library(ZOLTAN_LIBRARY zoltan
    PATHS "${ZOLTAN_PREFIX}/lib"
    NO_DEFAULT_PATH
  )
else()
  find_path(ZOLTAN_INCLUDE_DIR zoltan.h)
  find_library(ZOLTAN_LIBRARY zoltan)
endif()

set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY} )
set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR} )

if(ENABLE_PTSCOTCH)
  find_package(SCOTCH CONFIG REQUIRED)
elseif(ENABLE_PARMETIS)
  find_package(Parmetis MODULE REQUIRED)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    Zoltan
    DEFAULT_MSG
    ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR
)

mark_as_advanced(ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARY )
