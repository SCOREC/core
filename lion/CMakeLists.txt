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
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Do extra work if compression is enabled
if(LION_COMPRESS)
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(lion ${ZLIB_LIBRARIES})
endif()

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET lion
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS lion EXPORT lion-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT lion-targets DESTINATION lib/scorec)
