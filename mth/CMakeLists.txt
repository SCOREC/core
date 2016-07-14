# Package sources
set(MTH_SOURCES
  mthQR.cc
)

# Package headers
set(MTH_HEADERS
  mth.h
  mth_def.h
  mthVector.h
  mthMatrix.h
  mthTensor.h
  mthQR.h
  mthAD.h
)

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add the mth library
add_library(mth ${MTH_SOURCES})

# Link this library to these libraries
target_link_libraries(mth can)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET mth
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS mth EXPORT mth-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT mth-targets DESTINATION lib/scorec)
