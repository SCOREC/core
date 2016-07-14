# Package options
option(ENABLE_VIZ "Enable the viz package [ON|OFF]" OFF)
message(STATUS "ENABLE_VIZ: " ${ENABLE_VIZ})

# Only install the package if enabled
if(NOT ENABLE_VIZ)
  return()
endif()

# Find the milo library
find_path(MILO_INCLUDE_DIR NAMES milo.h)
if(NOT EXISTS "${MILO_INCLUDE_DIR}")
  message(FATAL_ERROR "milo header not found")
endif()
find_library(MILO_LIB NAMES milo)
if(NOT EXISTS "${MILO_LIB}")
 message(FATAL_ERROR "milo library not found")
endif()

# Package sources
set(SOURCES
  viz.cc
)

# Package headers
set(HEADERS
  viz.h
)

# Add the viz library
add_library(viz ${SOURCES})

# Inlcude directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MILO_INCLUDE_DIR}
)

# Link this package to these libraries
target_link_libraries(viz apf ${MILO_LIB})

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET viz
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS viz EXPORT viz-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT viz-targets DESTINATION lib/scorec)
