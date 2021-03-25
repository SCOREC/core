if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_EGADS)
  return()
endif()

# this file effectively sets the runtime path for the EGADS libraries
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gmi_egads_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/gmi_egads_config.h")

#Sources & Headers
set(SOURCES gmi_egads.c)
set(HEADERS gmi_egads.h)

add_library(gmi_egads ${SOURCES})

target_include_directories(gmi_egads
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${EGADS_INCLUDE_DIR}
)

target_link_libraries(gmi_egads
  PUBLIC
    gmi
    pcu
)

# Include directories
target_include_directories(gmi_egads PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
# make sure the compiler can find the config header
target_include_directories(gmi_egads PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

scorec_export_library(gmi_egads)

bob_end_subdir()
