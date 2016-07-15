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

# Add the mth library
add_library(mth ${MTH_SOURCES})

# Include directories
target_include_directories(mth INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this library to these libraries
# can is PUBLIC because public mth headers include can headers
target_link_libraries(mth PUBLIC can)

# Tell CMake to add this package to the SCOREC project
install(TARGETS mth EXPORT mth-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT mth-targets DESTINATION lib/scorec)
