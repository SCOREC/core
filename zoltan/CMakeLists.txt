if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(
  PUMI_ENABLE_ZOLTAN
  "Enable Zoltan interface [ON|OFF]"
  ${ENABLE_ZOLTAN} # Inherit from parent project
)
xsdk_add_tpl(ZOLTAN)
xsdk_add_tpl(PARMETIS)
message(STATUS "PUMI_ENABLE_ZOLTAN: " ${PUMI_ENABLE_ZOLTAN})

if(PUMI_NO_MPI AND PUMI_ENABLE_ZOLTAN)
  message(FATAL_ERROR "PUMI_NO_MPI is incompatible with Zoltan.")
endif()

# Let CMake know where to find custom find package implementations
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

# Check for and enable zoltan support
if(PUMI_ENABLE_ZOLTAN)
  if(Trilinos_PREFIX)
    set(ZOLTAN_PREFIX_DEFAULT "${Trilinos_PREFIX}")
    message(STATUS "ZOLTAN_PREFIX_DEFAULT ${ZOLTAN_PREFIX_DEFAULT}")
  endif()
  find_package(Zoltan MODULE REQUIRED)
endif()

# Package sources
if(PUMI_ENABLE_ZOLTAN)
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

# Package headers
set(HEADERS apfZoltan.h)

# Add the apf_zoltan library
add_library(apf_zoltan ${SOURCES})

# Include directories
target_include_directories(apf_zoltan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

# Link this package to these libraries
target_link_libraries(apf_zoltan PUBLIC pcu apf)

# Do extra work if zoltan is enabled
if(PUMI_ENABLE_ZOLTAN)
  target_include_directories(apf_zoltan PRIVATE
      ${ZOLTAN_INCLUDE_DIRS}
      ${PARMETIS_INCLUDE_DIRS}
      )
  target_link_libraries(apf_zoltan PUBLIC
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
  )
  target_compile_definitions(apf_zoltan PUBLIC PUMI_HAS_ZOLTAN)
endif()

scorec_export_library(apf_zoltan)

bob_end_subdir()
