# Package options

# Only install the package if enabled
if (NOT ENABLE_OMEGA_H)
  return()
endif()

# Package sources
set(SOURCES
  apfOmega_h.cc
)

# Package headers
set(HEADERS
  apfOmega_h.h
)

# Add the apf_omega_h library
add_library(apf_omega_h ${SOURCES})
target_compile_options(apf_omega_h PUBLIC "-std=c++11")

# Include directories
target_include_directories(apf_omega_h INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(apf_omega_h PUBLIC apf mds Omega_h::omega_h)

scorec_export_library(apf_omega_h)

bob_end_subdir()
