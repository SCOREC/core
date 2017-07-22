tribits_package(SCORECapf)

set(APF_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../mth
    ${CMAKE_CURRENT_SOURCE_DIR}/../can
    ${CMAKE_CURRENT_SOURCE_DIR}/../lion)

#Sources & Headers
set(APF_SOURCES
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
  apfMixedNumbering.cc
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
  apfSimplexAngleCalcs.cc
  apfFile.cc
)

set(APF_HEADERS
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
)

set(APF_SOURCES
    ${APF_SOURCES}
    ../mth/mthQR.cc
    ../lion/lionBase64.cc
    ../lion/lionNoZLib.cc
   )
set(APF_HEADERS
    ${APF_HEADERS}
    ../can/canArray.h
    ../can/canNewArray.h
    ../mth/mth.h
    ../mth/mth_def.h
    ../mth/mthVector.h
    ../mth/mthMatrix.h
    ../mth/mthTensor.h
    ../mth/mthQR.h
    ../mth/mthAD.h
   )

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${APF_INCLUDE_DIRS})

#Library
tribits_add_library(
   apf
   HEADERS ${APF_HEADERS}
   SOURCES ${APF_SOURCES})

tribits_package_postprocess()
