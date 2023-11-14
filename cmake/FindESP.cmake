#[=======================================================================[.rst:
FindESP
---------

Find ESP include dirs and libraries

Use this module by invoking :command:`find_package` with the form:

.. code-block:: cmake

  find_package(ESP
    [version] [EXACT]      # Minimum or EXACT version e.g. 1.19.0
    [REQUIRED]             # Fail with error if ESP is not found
    [COMPONENTS <libs>...] # ESP libraries by their canonical name
                           # e.g. "egads" for "libegads"
    [OPTIONAL_COMPONENTS <libs>...]
                           # Optional ESP libraries by their canonical name
  )                        # e.g. "egads" for "libegads"

This module finds headers and requested component libraries from ESP

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``ESP_FOUND``
  True if headers and requested libraries were found.

``ESP_INCLUDE_DIRS``
  ESP include directories.

``ESP_LIBRARY_DIRS``
  Link directories for ESP libraries.

``ESP_LIBRARIES``
  ESP component libraries to be linked.

``ESP_<COMPONENT>_FOUND``
  True if component ``<COMPONENT>`` was found.

``ESP_<COMPONENT>_LIBRARY``
  Libraries to link for component ``<COMPONENT>`` (may include
  :command:`target_link_libraries` debug/optimized keywords).

Cache variables
^^^^^^^^^^^^^^^

Search results are saved persistently in CMake cache entries:

``ESP_INCLUDE_DIR``
  Directory containing ESP headers.

``ESP_LIBRARY_DIR``
  Directory containing ESP libraries.

Hints
^^^^^

This module reads hints about search locations from variables:

``ESP_ROOT``, ``ESPROOT``
  Preferred installation prefix.

``ESP_INCLUDEDIR``
  Preferred include directory e.g. ``<prefix>/include``.

``ESP_LIBRARYDIR``
  Preferred library directory e.g. ``<prefix>/lib``.

``ESP_NO_SYSTEM_PATHS``
  Set to ``ON`` to disable searching in locations not
  specified by these hint variables. Default is ``OFF``.

``ESP_ADDITIONAL_VERSIONS``
  List of ESP versions not known to this module.
  (ESP install locations may contain the version).

Users may set these hints or results as ``CACHE`` entries.  Projects
should not read these entries directly but instead use the above
result variables.  Note that some hint names start in upper-case
``ESP``.  One may specify these as environment variables if they are
not specified as CMake variables or cache entries.

This module first searches for the ESP header files using the above
hint variables (excluding ``ESP_LIBRARYDIR``) and saves the result in
``ESP_INCLUDE_DIR``.  Then it searches for requested component libraries
using the above hints (excluding ``ESP_INCLUDEDIR``), "lib" directories 
near ``ESP_INCLUDE_DIR``, and the library name configuration settings below. 
It saves the library directories in ``ESP_LIBRARY_DIR`` and individual library
locations in ``ESP_<COMPONENT>_LIBRARY``.
When one changes settings used by previous searches in the same build
tree (excluding environment variables) this module discards previous
search results affected by the changes and searches again.

Imported Targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``ESP::ESP``
  Interface target for all components linking against all components.

``ESP::<component>``
  Target for specific component dependency (shared or static library).

It is important to note that the imported targets behave differently
than variables created by this module: multiple calls to
:command:`find_package(ESP)` in the same directory or sub-directories with
different options (e.g. static or shared) will not override the
values of the targets created by the first call.

Examples
^^^^^^^^

Find ESP libraries and use imported targets:

.. code-block:: cmake

  find_package(ESP REQUIRED COMPONENTS
               egads ocsm)
  add_executable(foo foo.cc)
  target_link_libraries(foo ESP::egads ESP::aimUtil)

#]=======================================================================]

include(GNUInstallDirs)

set(quiet "")
if(ESP_FIND_QUIETLY)
  set(quiet QUIET)
endif()

# ------------------------------------------------------------------------
# Find ESP include dir
# ------------------------------------------------------------------------
if(NOT ESP_INCLUDE_DIR)

  set(_ESP_INCLUDE_SEARCH_DIRS "")
  if(ESP_INCLUDEDIR)
    list(APPEND _ESP_INCLUDE_SEARCH_DIRS ${ESP_INCLUDEDIR})
  endif()

  if(DEFINED ENV{ESP_ROOT})
    list(APPEND _ESP_INCLUDE_SEARCH_DIRS $ENV{ESP_ROOT}/include $ENV{ESP_ROOT})
  endif()

  if(DEFINED ENV{ESPROOT})
    list(APPEND _ESP_INCLUDE_SEARCH_DIRS $ENV{ESPROOT}/include $ENV{ESPROOT})
  endif()

  find_path(ESP_INCLUDE_DIR NAMES egads.h HINTS ${_ESP_INCLUDE_SEARCH_DIRS})
endif()

message(STATUS "ESP include dir: ${ESP_INCLUDE_DIR}")

# ------------------------------------------------------------------------
#  Extract version information from egadsTypes.h
# ------------------------------------------------------------------------
if(ESP_INCLUDE_DIR)

  # Extract ESP_VERSION_MAJOR AND ESP_VERISON_MINOR from egadsTypes.h
  set(ESP_VERSION_MAJOR 0)
  set(ESP_VERSION_MINOR 0)
  file(STRINGS "${ESP_INCLUDE_DIR}/egadsTypes.h" _ESP_VERSION_CONTENTS REGEX "#define EGADSMAJOR ")
  if("${_ESP_VERSION_CONTENTS}" MATCHES "#define EGADSMAJOR[ \t\r\n]+([0-9]+)")
    set(ESP_VERSION_MAJOR "${CMAKE_MATCH_1}")
  endif()
  unset(_ESP_VERSION_HEADER_CONTENTS)

  file(STRINGS "${ESP_INCLUDE_DIR}/egadsTypes.h" _ESP_VERSION_CONTENTS REGEX "#define EGADSMINOR ")
  if("${_ESP_VERSION_CONTENTS}" MATCHES "#define EGADSMINOR[ \t\r\n]+([0-9]+)")
    set(ESP_VERSION_MINOR "${CMAKE_MATCH_1}")
  endif()
  unset(_ESP_VERSION_HEADER_CONTENTS)

  # ESP versioning does not include a patch number so we set it to zero
  SET(ESP_VERSION_PATCH 0)

  # Define alias variables for backwards compat.
  set(ESP_MAJOR_VERSION ${ESP_VERSION_MAJOR})
  set(ESP_MINOR_VERSION ${ESP_VERSION_MINOR})
  set(ESP_SUBMINOR_VERSION ${ESP_VERSION_PATCH})

  # Define ESP version in x.y.z format
  set(ESP_VERSION_STRING "${ESP_VERSION_MAJOR}.${ESP_VERSION_MINOR}.${ESP_VERSION_PATCH}")

  # message(STATUS "ESP Version: ${ESP_VERSION_STRING}")
endif()

# ------------------------------------------------------------------------
#  Begin finding ESP libraries
# ------------------------------------------------------------------------

# all potential ESP components
set(ESP_COMPONENTS caps egads ocsm)

# if not explicitly asking for any component, find all of them
if(NOT ESP_FIND_COMPONENTS)
  set(ESP_FIND_COMPONENTS ${ESP_COMPONENTS})
endif()

foreach(component ${ESP_FIND_COMPONENTS})
  
  if(component STREQUAL "egads")
    if (ESP_USE_EGADSLITE)
      find_library(egads_LIBRARY NAMES egadslite)
    else()
      find_library(egads_LIBRARY NAMES egads)  
    endif()
  else()
    find_library(${component}_LIBRARY NAMES ${component})
  endif()

  if(${component}_LIBRARY)
    set(ESP_${component}_FOUND True)
  else()
    set(ESP_${component}_FOUND False)
  endif()

  # Create a library target only if the above checks passed
  if(ESP_${component}_FOUND AND NOT TARGET ESP::${component})
    # Can't easily tell how ESP was compiled, so just default to UNKNOWN
    # library type and CMake will make a best effort guess
    add_library(ESP::${component} UNKNOWN IMPORTED)

    set_property(
        TARGET ESP::${component} PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${ESP_INCLUDE_DIR}"
    )
    if(EXISTS "${${component}_LIBRARY}")
      set_property(
          TARGET ESP::${component} PROPERTY
          IMPORTED_LOCATION "${${component}_LIBRARY}"
      )
    endif()
  endif()
endforeach()

# Create INTERFACE target that bundles all the found libraries together
if(NOT TARGET ESP::ESP)
  add_library(ESP::ESP INTERFACE IMPORTED)
  foreach(component ${ESP_FIND_COMPONENTS})
    if(TARGET ESP::${component})
      target_link_libraries(ESP::ESP INTERFACE ESP::${component})
    endif()
  endforeach()
endif()

# Use CMake provided module to check the variables
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ESP
  REQUIRED_VARS ESP_INCLUDE_DIR
  VERSION_VAR ESP_VERSION_STRING
  HANDLE_COMPONENTS
)