if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_CAPSTONE)
  return()
endif()

#Sources & Headers
set(SOURCES gmi_cap.cc)
set(HEADERS gmi_cap.h)

add_library(gmi_cap ${SOURCES})

target_link_libraries(gmi_cap PUBLIC gmi pcu capstone_module framework_testing)

# Include directories
target_include_directories(gmi_cap PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
# make sure the compiler can find the config header
target_include_directories(gmi_cap PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

scorec_export_library(gmi_cap)

bob_end_subdir()
