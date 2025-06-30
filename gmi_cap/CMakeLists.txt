if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT PUMI_ENABLE_CAPSTONE)
  return()
endif()

#Sources & Headers
set(SOURCES gmi_cap.cc)
set(HEADERS gmi_cap.h)

add_library(gmi_cap ${SOURCES})

target_link_libraries(gmi_cap PUBLIC gmi framework_core)
target_link_libraries(gmi_cap PRIVATE pcu lion capstone_module)

# Include directories
target_include_directories(gmi_cap PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )
# make sure the compiler can find the config header
target_include_directories(gmi_cap PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

scorec_export_library(gmi_cap)

bob_end_subdir()
