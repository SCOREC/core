# Package options
option(ENABLE_DSP "Enable the dsp package [ON|OFF]" OFF)
message(STATUS "ENABLE_DSP: " ${ENABLE_DSP})

# Only install the package if enabled
if(NOT ENABLE_DSP)
  return()
endif()

# Package sources
set(SOURCES
  dspAdapters.cc
  dspSmoothers.cc
  dspGraphDistance.cc
  dsp.cc
)

# Package headers
set(HEADERS
  dspAdapters.h
  dspSmoothers.h
  dsp.h
)

# Add the dsp library
add_library(dsp ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link this package to these libraries
target_link_libraries(dsp ma)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET dsp
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS dsp EXPORT dsp-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT dsp-targets DESTINATION lib/scorec)
