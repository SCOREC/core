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
  crvShapeFixer.cc
  crvShapeHandler.cc
  crvSnap.cc
  crvSubdivision.cc
  crvTables.cc
  crvQuality.cc
  crvVtk.cc
  LBFGS.cc
  crvEdgeOptim.cc
  crvModelEdgeOptim.cc
  crvFaceOptim.cc
  optim.cc
)

# Package headers
set(HEADERS
  crv.h
  crvEdgeOptim.h
  crvModelEdgeOptim.h
  crvFaceOptim.h
  LBFGS.h
  optim.h
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
