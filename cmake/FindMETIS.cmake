# Copyright (C) 2025 Scientific Computation Research Center
#
# This work is open source software, licensed under the terms of the
# BSD license as described in the LICENSE file in the top-level directory.
#
#[====================================[
FindMETIS.cmake

Hints:
- `METIS_PREFIX`: The prefix for the METIS installation.

Once done this will define:
- `METIS_FOUND`: True if METIS is found.
- `METIS::METIS`: The IMPORTED library if found.
- `METIS_LIBRARIES`: METIS libraries required for linking.
- `METIS_INCLUDE_DIRS`: Directories containing METIS headers.
]====================================]

cmake_policy(PUSH)

set(METIS_PREFIX "${Trilinos_PREFIX}" CACHE STRING "METIS install directory")
find_path(METIS_INCLUDE_DIR metis.h HINTS "${METIS_PREFIX}/include")
set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _METIS_VERSION_STRS
  REGEX "#define METIS_VER_(MAJOR|(SUB)?MINOR)"
)
set(METIS_VERSION_STR "")
foreach(ver_comp IN LISTS _METIS_VERSION_STRS)
  if(METIS_VERSION_STR)
    string(APPEND METIS_VERSION_STR .)
  endif()
  string(REGEX REPLACE
    "^#define METIS_VER_(MAJOR|(SUB)?MINOR)[ \t]+" ""
    comp_num "${ver_comp}"
  )
  string(APPEND METIS_VERSION_STR "${comp_num}")
endforeach()
find_library(METIS_LIBRARY metis HINTS "${METIS_PREFIX}/lib")
# Add imported library.
add_library(METIS::METIS UNKNOWN IMPORTED)
set_target_properties(METIS::METIS PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}"
  IMPORTED_LINK_INTERFACE_LANGUAGE "C;CXX"
  IMPORTED_LOCATION "${METIS_LIBRARY}"
)
# Sometimes GKLib is an external dependency; usually it is in METIS itself.
find_library(GK_LIBRARY GKlib HINTS "${METIS_PREFIX}/lib")
if(EXISTS "${GK_LIBRARY}")
  set_target_property(METIS::METIS PROPERTIES
    INTERFACE_LINK_LIBRARIES "${GK_LIBRARY}"
  )
  set(METIS_LIBRARIES "${METIS_LIBRARY}" "${GK_LIBRARY}")
else()
  set(METIS_LIBRARIES "${METIS_LIBRARY}")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY
  VERSION_VAR METIS_VERSION_STR
)
mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY)
cmake_policy(POP)
