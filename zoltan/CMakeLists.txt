if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(ENABLE_ZOLTAN "Enable Zoltan interface [ON|OFF]" OFF)
xsdk_add_tpl(ZOLTAN)
xsdk_add_tpl(PARMETIS)
message(STATUS "ENABLE_ZOLTAN: " ${ENABLE_ZOLTAN})

# Let CMake know where to find custom find package implementations
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

# Check for and enable zoltan support
if(ENABLE_ZOLTAN)
  if(Trilinos_PREFIX)
    set(ZOLTAN_PREFIX_DEFAULT "${Trilinos_PREFIX}")
    message(STATUS "ZOLTAN_PREFIX_DEFAULT ${ZOLTAN_PREFIX_DEFAULT}")
  endif()
  find_package(Zoltan MODULE REQUIRED)
endif()

# Package sources
if(ENABLE_ZOLTAN)
  set(SOURCES
    apfInterElement.cc
    apfZoltan.cc
    apfZoltanMesh.cc
    apfZoltanCallbacks.cc
  )
else()
  set(SOURCES
    apfInterElement.cc
    apfZoltanEmpty.cc
  )
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/apfZoltanConfig.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/apfZoltanConfig.h")

# Set include directory for building apf_zoltan target
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# Package headers
set(HEADERS apfZoltan.h "${CMAKE_CURRENT_BINARY_DIR}/apfZoltanConfig.h")

# Add the apf_zoltan library
add_library(apf_zoltan ${SOURCES})

# Include directories
target_include_directories(apf_zoltan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Set include directory for generated header for targets linking to apf_zoltan
target_include_directories(apf_zoltan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    )


# Link this package to these libraries
target_link_libraries(apf_zoltan PUBLIC pcu apf)

# Do extra work if zoltan is enabled
if(ENABLE_ZOLTAN)
  target_include_directories(apf_zoltan PRIVATE
      ${ZOLTAN_INCLUDE_DIRS}
      ${PARMETIS_INCLUDE_DIRS}
      )
  target_link_libraries(apf_zoltan PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
  )
endif()

scorec_export_library(apf_zoltan)

bob_end_subdir()
