# Component options
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
target_include_directories(apf INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

# Link this package to these libraries
target_link_libraries(dsp PUBLIC ma apf)

scorec_export_library(dsp)

bob_end_subdir()
