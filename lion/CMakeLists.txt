# Package options
option(LION_COMPRESS "Enable zlib compression [ON|OFF]" OFF)
message(STATUS "LION_COMPRESS: " ${LION_COMPRESS})

# Check for and enable zlib support
if (LION_COMPRESS)
  find_package(ZLIB REQUIRED)
endif()

# Package sources
set(SOURCES lionBase64.cc)
if(LION_COMPRESS)
  set(SOURCES ${SOURCES} lionZLib.cc)
else()
  set(SOURCES ${SOURCES} lionNoZLib.cc)
endif()

# Package headers
set(HEADERS
  lionBase64.h
  lionCompress.h
)

# Add the lion Library
add_library(lion ${SOURCES})

# Include directories
target_include_directories(lion INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Do extra work if compression is enabled
if(LION_COMPRESS)
  target_include_directories(lion PRIVATE ${ZLIB_INCLUDE_DIR})
  target_link_libraries(lion PUBLIC ${ZLIB_LIBRARIES})
endif()

bob_export_target(lion)

bob_end_subdir()
