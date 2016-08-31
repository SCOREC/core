# Package sources
set(SOURCES
  mthQR.cc
)

# Package headers
set(HEADERS
  mth.h
  mth_def.h
  mthVector.h
  mthMatrix.h
  mthTensor.h
  mthQR.h
  mthAD.h
)

# Add the mth library
add_library(mth ${SOURCES})

# Include directories
target_include_directories(mth INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this library to these libraries
target_link_libraries(mth PUBLIC can)

scorec_export_library(mth)

bob_end_subdir()
