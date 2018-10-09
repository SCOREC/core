tribits_package(SCOREClion)

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

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Library
tribits_add_library(
   lion
   HEADERS ${HEADERS}
   SOURCES ${SOURCES})

# Do extra work if compression is enabled
if(LION_COMPRESS)
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(lion ${ZLIB_LIBRARIES})
endif()

tribits_package_postprocess()
