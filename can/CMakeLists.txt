# Package headers
set(HEADERS
  canArray.h
  canNewArray.h
)

# Add the can 'library' as a CMake Interface
# (header only library)
add_library(can INTERFACE)

# Include directories
target_include_directories(can INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

scorec_export_library(can)

bob_end_subdir()
