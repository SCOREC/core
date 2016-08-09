# Package sources
set(SOURCES
  sprGetGradIPField.cc
  sprRecoverField.cc
  sprEstimateError.cc
  sprEstimateTargetError.cc)

# Package headers
set(HEADERS
  spr.h)

# Add the spr library
add_library(spr ${SOURCES})

# Include directories
target_include_directories(spr INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(spr PUBLIC apf PRIVATE pcu)

target_macro(spr)
