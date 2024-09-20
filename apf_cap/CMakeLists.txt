if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_CAPSTONE)
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

if(HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_compile_definitions(apf_cap PRIVATE HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_link_libraries(apf_cap PRIVATE framework_meshing)
target_compile_features(apf_cap PRIVATE cxx_std_14)
endif()

include(GNUInstallDirs)

target_include_directories(apf_cap PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

scorec_export_library(apf_cap)

bob_end_subdir()
