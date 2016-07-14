# Package options
option(ENABLE_VIZ "Enable the viz package [ON|OFF]" OFF)
message(STATUS "ENABLE_VIZ: " ${ENABLE_VIZ})

# Only install the package if enabled
if(NOT ENABLE_VIZ)
  return()
endif()

# Find the milo library
find_path(MILO_INCLUDE NAMES milo.h)
find_library(MILO_LIB NAMES milo)

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
  ${MILO_INCLUDE}
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
