set(SOURCES
  samSz.cc
  sam.cc
  samElementCount.cc
)

# Package headers
set(HEADERS
  sam.h
  samSz.h
  samElementCount.h
)

# Add the sam library
add_library(sam ${SOURCES})

# Include directories
target_include_directories(sam INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(sam
    PUBLIC apf
    PRIVATE pcu)

target_macro(sam)
