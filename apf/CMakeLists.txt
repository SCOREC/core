if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package sources
set(SOURCES
  apf.cc
  apfCavityOp.cc
  apfElement.cc
  apfField.cc
  apfFieldOf.cc
  apfGradientByVolume.cc
  apfIntegrate.cc
  apfMatrix.cc
  apfDynamicMatrix.cc
  apfDynamicVector.cc
  apfMatrixField.cc
  apfMatrixElement.cc
  apfMesh.cc
  apfMesh2.cc
  apfMigrate.cc
  apfScalarElement.cc
  apfScalarField.cc
  apfShape.cc
  apfIPShape.cc
  apfHierarchic.cc
  apfPolyBasis1D.cc
  apfNedelec.cc
  apfL2Shapes.cc
  apfH1Shapes.cc
  apfVector.cc
  apfVectorElement.cc
  apfVectorField.cc
  apfMixedVectorElement.cc
  apfMixedVectorField.cc
  apfPackedField.cc
  apfNumbering.cc
  apfMixedNumbering.cc
  apfAdjReorder.cc
  apfVtk.cc
  apfVtkPieceWiseFields.cc
  apfFieldData.cc
  apfTagData.cc
  apfCoordData.cc
  apfArrayData.cc
  apfUserData.cc
  apfPartition.cc
  apfConvert.cc
  apfConstruct.cc
  apfVerify.cc
  apfGeometry.cc
  apfBoundaryToElementXi.cc
  apfSimplexAngleCalcs.cc
  apfFile.cc
  apfMIS.cc
)

if(ENABLE_CGNS)
  set(SOURCES ${SOURCES} apfCGNS.cc)
endif(ENABLE_CGNS)

# Package headers
set(HEADERS
  apf.h
  apfMesh.h
  apfMesh2.h
  apfMatrix.h
  apfVector.h
  apfArray.h
  apfDynamicMatrix.h
  apfDynamicVector.h
  apfDynamicArray.h
  apfNew.h
  apfCavityOp.h
  apfShape.h
  apfNumbering.h
  apfMixedNumbering.h
  apfPartition.h
  apfConvert.h
  apfGeometry.h
  apf2mth.h
  apfMIS.h
  apfField.h
  apfFieldData.h
  apfNumberingClass.h
  apfElement.h
)

# Add the apf library
add_library(apf ${SOURCES})

# Include directories
target_include_directories(apf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this library to these others
target_link_libraries(apf
   PUBLIC
     pcu
     gmi
     lion
     can
     mth
   )

scorec_export_library(apf)

bob_end_subdir()
