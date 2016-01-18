cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

project(apf_omega_h)

option(ENABLE_OMEGA_H "Enable Omega_h interface" OFF)

if (NOT ENABLE_OMEGA_H)
  return()
endif()

find_package(mds PATHS ${CMAKE_BINARY_DIR})
find_package(apf PATHS ${CMAKE_BINARY_DIR})

find_path(OMEGA_H_INCLUDE_DIR omega_h.h)
if(NOT EXISTS "${OMEGA_H_INCLUDE_DIR}")
  message(FATAL_ERROR "Omega_h header not found")
endif()
find_library(OMEGA_H_LIB omega_h)
if(NOT EXISTS "${OMEGA_H_LIB}")
 message(FATAL_ERROR "Omega_h library not found")
endif()

set(DEP_LIBS ${MDS_LIBS} ${APF_LIBS} ${OMEGA_H_LIB})
set(APF_OMEGA_H_LIBS apf_omega_h ${DEP_LIBS})
set(APF_OMEGA_H_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MDS_INCLUDE_DIRS}
  ${APF_INCLUDE_DIRS}
  ${OMEGA_H_INCLUDE_DIR})

set(APF_OMEGA_H_SOURCES apfOmega_h.cc)
set(APF_OMEGA_H_HEADERS apfOmega_h.h)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${APF_OMEGA_H_INCLUDE_DIRS})

#Library
add_library(apf_omega_h ${APF_OMEGA_H_SOURCES})
target_link_libraries(apf_omega_h ${DEP_LIBS})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_omega_hConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apf_omega_hConfig.cmake")

INSTALL(FILES ${APF_OMEGA_H_HEADERS} DESTINATION include)
INSTALL(TARGETS
  apf_omega_h
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
