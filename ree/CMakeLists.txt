if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package sources
set(SOURCES
  reeResidualFunctionals.cc
  reeFluxCorrection.cc
  reeCorrectedFlux.cc
  reeEstimateError.cc
  reeSizeField.cc)

# Package headers
set(HEADERS
  ree.h)

# Add the ree library
add_library(ree ${SOURCES})

# Include directories
target_include_directories(ree INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(ree PUBLIC apf pcu crv)

scorec_export_library(ree)

bob_end_subdir()
