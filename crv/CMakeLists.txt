# Package sources
set(SOURCES
  crv.cc
  crvDBG.cc
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
  crvShapeFixer.cc
  crvShapeHandler.cc
  crvSnap.cc
  crvSubdivision.cc
  crvTables.cc
  crvQuality.cc
  crvVtk.cc
  crvOptimizations.cc
  crvObjectiveFunctions.cc
)

# Package headers
set(HEADERS
  crv.h
  crvDBG.h
  crvOptimizations.h
  crvObjectiveFunctions.h
  crvBezier.h
  crvTables.h
  crvQuality.h
)

# Add the crv library
add_library(crv ${SOURCES})

# Include directories
target_include_directories(crv INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# TODO: Add these to the higher lever cmake file in a proper form
target_include_directories(crv PUBLIC
  "/lore/hakimm2/opt/eigen/include/eigen3"
  )
target_include_directories(crv PUBLIC
  "/lore/hakimm2/opt/LBFGSpp/include"
  )
# Link this package to these libraries
target_link_libraries(crv
    PUBLIC
      apf
      ma
      gmi
      pcu
    )

scorec_export_library(crv)

bob_end_subdir()
