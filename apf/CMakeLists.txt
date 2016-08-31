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
  apfMesh.cc
  apfMesh2.cc
  apfMigrate.cc
  apfScalarElement.cc
  apfScalarField.cc
  apfShape.cc
  apfIPShape.cc
  apfHierarchic.cc
  apfVector.cc
  apfVectorElement.cc
  apfVectorField.cc
  apfPackedField.cc
  apfNumbering.cc
  apfAdjReorder.cc
  apfVtk.cc
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
  apfFile.cc
)

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
  apfPartition.h
  apfConvert.h
  apfGeometry.h
  apf2mth.h
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
