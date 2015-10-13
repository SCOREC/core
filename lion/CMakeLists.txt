cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

project(lion)

set(LION_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR})

set(LION_LIBS lion)

set(LION_SOURCES
)

set(LION_HEADERS
)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${LION_INCLUDE_DIRS})

#Library
add_library(lion ${LION_SOURCES})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lionConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/lionConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/liblion.pc.in"
    "${PROJECT_BINARY_DIR}/liblion.pc"
    @ONLY)
#pkgconfig }

INSTALL(FILES ${MTH_HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/liblion.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS
  lion
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
