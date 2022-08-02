if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_EGADS)
  return()
endif()

#Sources & Headers
set(SOURCES gmi_egads.c)
set(HEADERS gmi_egads.h)

add_library(gmi_egads ${SOURCES})

target_link_libraries(gmi_egads
  PUBLIC
    gmi
  PRIVATE
    ESP::egads
)

target_compile_definitions(gmi_egads
  INTERFACE
    PUMI_HAS_EGADS
)

if (PUMI_USE_EGADSLITE)
  target_compile_definitions(gmi_egads
    PRIVATE
      PUMI_USE_EGADSLITE
  )
endif()

# Include directories
target_include_directories(gmi_egads PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

scorec_export_library(gmi_egads)

bob_end_subdir()
