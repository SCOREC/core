if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package sources
set(SOURCES
  gmi.c
  agm.c
  gmi_base.c
  gmi_file.c
  gmi_lookup.c
  gmi_mesh.c
  gmi_null.c
  gmi_analytic.c
)

# Package headers
set(HEADERS
  gmi.h
  agm.h
  gmi_base.h
  gmi_lookup.h
  gmi_mesh.h
  gmi_null.h
  gmi_analytic.h
)

# Add the gmi library
add_library(gmi ${SOURCES})
target_link_libraries(gmi pcu lion)
target_include_directories(gmi INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

scorec_export_library(gmi)

bob_end_subdir()
