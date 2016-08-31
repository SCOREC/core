if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(ENABLE_ZOLTAN "Enable Zoltan interface [ON|OFF]" OFF)
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

# Package headers
set(HEADERS apfZoltan.h)

# Add the apf_zoltan library
add_library(apf_zoltan ${SOURCES})

# Include directories
target_include_directories(apf_zoltan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
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
