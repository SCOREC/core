cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "0.1")

project(viz)

option(ENABLE_VIZ "Enable the viz package [ON|OFF]" OFF)
if (ENABLE_VIZ)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_path(MILO_INCLUDE NAMES milo.h)
  find_library(MILO_LIB NAMES milo)
  include_directories(${APF_INCLUDE_DIRS})
  include_directories(${MILO_INCLUDE})
  set(VIZ_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
  set(DEP_LIBS ${APF_LIBS} ${MILO_LIB})
  set(VIZ_LIBS viz ${DEP_LIBS})
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/vizConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/vizConfig.cmake")

if (NOT ENABLE_VIZ)
  return()
endif()

include_directories(${VIZ_INCLUDE_DIRS})

#Sources & Headers
set(SOURCES
   viz.cc
)

set(HEADERS
   viz.h
)

#Library
add_library(viz ${SOURCES})
target_link_libraries(viz ${DEP_LIBS})

#Install
INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(TARGETS  
   viz
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)

