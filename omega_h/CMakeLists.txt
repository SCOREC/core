# Package options
option(ENABLE_OMEGA_H "Enable the Omega_h interface" OFF)
message(STATUS "ENABLE_OMEGA_H:" ${ENABLE_OMEGA_H})

# Only install the package if enabled
if (NOT ENABLE_OMEGA_H)
  return()
endif()

# find the omega_h library
find_path(OMEGA_H_INCLUDE_DIR omega_h.h)
if(NOT EXISTS "${OMEGA_H_INCLUDE_DIR}")
  message(FATAL_ERROR "Omega_h header not found")
endif()
find_library(OMEGA_H_LIB omega_h)
if(NOT EXISTS "${OMEGA_H_LIB}")
 message(FATAL_ERROR "Omega_h library not found")
endif()

# Package sources
set(SOURCES
  apfOmega_h.cc
)

# Package headers
set(SOURCES
  apfOmega_h.h
)

# Add the apf_omega_h library
add_library(apf_omega_h ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${OMEGA_H_INCLUDE_DIR}
)

# Link this package to these libraries
target_link_libraries(apf_omega_h mds apf ${OMEGA_H_LIB})

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET apf_omega_h
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS apf_omega_h EXPORT apf_omega_h-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT apf_omega_h-targets DESTINATION lib/scorec)
