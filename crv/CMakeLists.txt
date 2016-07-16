# Package sources
set(SOURCES
  crv.cc
  crvAdapt.cc
  crvBernstein.cc
  crvBezier.cc
  crvBezierPoints.cc
  crvBezierShapes.cc
  crvBlended.cc
  crvCurveMesh.cc
  crvElevation.cc
  crvG1Points.cc
  crvMath.cc
  crvReposition.cc
  crvShape.cc
  crvShapeHandler.cc
  crvSnap.cc
  crvSubdivision.cc
  crvTables.cc
  crvQuality.cc
  crvVtk.cc
)

# Package headers
set(HEADERS
  crv.h
)

# Add the crv library
add_library(crv ${SOURCES})

# Include directories
target_include_directories(ma INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(crv PUBLIC apf ma)

# Tell CMake to add this package to the SCOREC project
install(TARGETS crv EXPORT crv-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT crv-targets DESTINATION lib/scorec)
