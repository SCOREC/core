cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

project(lion)

option(LION_COMPRESS "Enable zlib compression [ON|OFF]" OFF)
message(STATUS "LION_COMPRESS " ${LION_COMPRESS})
if (LION_COMPRESS)
  find_package(ZLIB REQUIRED)
  message(STATUS "ZLIB_FOUND " ${ZLIB_FOUND})
  message(STATUS "ZLIB_INCLUDE_DIR " ${ZLIB_INCLUDE_DIR})
  message(STATUS "ZLIB_LIBRARIES " ${ZLIB_LIBRARIES})
  include_directories(${ZLIB_INCLUDE_DIR})
  set(DEP_LIBS ${DEP_LIBS} ${ZLIB_LIBRARIES})
  message(STATUS "DEP_LIBS " ${DEP_LIBS})
endif()

set(LION_LIBS lion ${DEP_LIBS})
set(LION_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR})

set(LION_SOURCES lionBase64.cc)
set(LION_HEADERS lionBase64.h lionCompress.h)

if (LION_COMPRESS)
  set(LION_SOURCES ${LION_SOURCES} lionZLib.cc)
else()
  set(LION_SOURCES ${LION_SOURCES} lionNoZLib.cc)
endif()

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${LION_INCLUDE_DIRS})

#Library
add_library(lion ${LION_SOURCES})
target_link_libraries(lion ${DEP_LIBS})
message(STATUS "target_link_libraries " ${DEP_LIBS})

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

INSTALL(FILES ${LION_HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/liblion.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS
  lion
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
