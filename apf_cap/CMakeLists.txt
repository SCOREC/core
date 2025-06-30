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
set(SOURCES apfCAP.cc apfCAPsizing.cc)
set(HEADERS apfCAP.h)

#Library
add_library(apf_cap ${SOURCES})
target_link_libraries(apf_cap PUBLIC apf gmi_cap)
target_link_libraries(apf_cap PRIVATE framework_mesh framework_application)

set(CMAKE_CXX_OLD_STANDARD "${CMAKE_CXX_STANDARD}")
cmake_push_check_state(RESET)
set(CMAKE_CXX_STANDARD 14)

# CreateFT_Vector.h via CreateMG_SizingMetricTool.h is missing
# `#include <limits>` in Capstone 12.5.2 so add it explicitly in this check.
if(MSVC)
  set(CMAKE_REQUIRED_FLAGS "/FI limits")
else()
  set(CMAKE_REQUIRED_FLAGS "-include limits")
endif()

set(CMAKE_REQUIRED_LIBRARIES "framework_meshing")
check_include_file_cxx("CreateMG_SizingMetricTool.h" HAVE_CAPSTONE_SIZINGMETRICTOOL)
cmake_pop_check_state()
set(CMAKE_CXX_STANDARD "${CMAKE_CXX_OLD_STANDARD}")

if(HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_compile_definitions(apf_cap PRIVATE HAVE_CAPSTONE_SIZINGMETRICTOOL)
target_compile_features(framework_meshing INTERFACE cxx_std_14)
target_link_libraries(apf_cap PRIVATE framework_meshing)
endif()

include(GNUInstallDirs)

target_include_directories(apf_cap PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

scorec_export_library(apf_cap)

bob_end_subdir()
