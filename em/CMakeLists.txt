if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package sources
set(SOURCES
  emResidualFunctionals.cc
  emFluxCorrection.cc
  emEstimateError.cc)

# Package headers
set(HEADERS
  em.h)

# Add the em library
add_library(em ${SOURCES})

# Include directories
target_include_directories(em INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(em PUBLIC apf pcu crv)

scorec_export_library(em)

bob_end_subdir()
