if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_CAPSTONE)
  return()
endif()

include(CMakePushCheckState)
include(CheckIncludeFileCXX)
cmake_policy(SET CMP0075 NEW) # Observe CMAKE_REQUIRED_LIBRARIES.

#Sources & Headers
set(SOURCES apfCAP.cc)
set(HEADERS apfCAP.h)

#Library
add_library(apf_cap ${SOURCES})
target_link_libraries(apf_cap PUBLIC apf gmi_cap)
target_link_libraries(apf_cap PUBLIC capstone_module
				     framework_testing)

set(CMAKE_CXX_OLD_STANDARD "${CMAKE_CXX_STANDARD}")
cmake_push_check_state(RESET)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_REQUIRED_LIBRARIES "framework_meshing")
check_include_file_cxx("CreateMG_SizingMetricTool.h" HAVE_CAPSTONE_SIZINGMETRICTOOL)
cmake_pop_check_state()
set(CMAKE_CXX_STANDARD "${CMAKE_CXX_OLD_STANDARD}")

if(HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_compile_definitions(apf_cap PRIVATE HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_link_libraries(apf_cap PRIVATE framework_meshing)
target_compile_features(apf_cap PRIVATE cxx_std_14)
endif()

include(GNUInstallDirs)

target_include_directories(apf_cap PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

scorec_export_library(apf_cap)

bob_end_subdir()
