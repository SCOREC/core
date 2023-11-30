if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

#Sources & Headers
set(SOURCES apfCAP.cc)
set(HEADERS apfCAP.h)

#Library
add_library(apf_cap ${SOURCES})
target_link_libraries(apf_cap PUBLIC apf gmi_cap)
target_link_libraries(apf_cap PUBLIC capstone_module
				     framework_testing)
target_include_directories(apf_cap PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

#directory containing apf_simConfig.h
target_include_directories(apf_cap PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

scorec_export_library(apf_cap)

bob_end_subdir()
