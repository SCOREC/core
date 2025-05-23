if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(LION_COMPRESS "Enable zlib compression [ON|OFF]" OFF)
message(STATUS "LION_COMPRESS: " ${LION_COMPRESS})

# Check for and enable zlib support
if (LION_COMPRESS)
  find_package(ZLIB REQUIRED)
endif()

# Package sources
set(SOURCES lionBase64.cc lionPrint.c)
if(LION_COMPRESS)
  set(SOURCES ${SOURCES} lionZLib.cc)
else()
  set(SOURCES ${SOURCES} lionNoZLib.cc)
endif()

# Package headers
set(HEADERS
  lionBase64.h
  lionCompress.h
  lionPrint.h
)

# Add the lion Library
add_library(lion ${SOURCES})

# Include directories
target_include_directories(lion INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

# Do extra work if compression is enabled
if(LION_COMPRESS)
  target_include_directories(lion PRIVATE ${ZLIB_INCLUDE_DIR})
  target_link_libraries(lion PUBLIC ${ZLIB_LIBRARIES})
endif()

scorec_export_library(lion)

bob_end_subdir()
