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
    )

# Link this package to these libraries
target_link_libraries(sam PUBLIC apf pcu)

scorec_export_library(sam)

bob_end_subdir()
