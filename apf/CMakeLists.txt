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

#quickfix hack which is now semi-permanent
if (BUILD_IN_TRILINOS)
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
endif()

# Add the apf library
add_library(apf ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu
  ${CMAKE_CURRENT_SOURCE_DIR}/../gmi
  ${CMAKE_CURRENT_SOURCE_DIR}/../can
  ${CMAKE_CURRENT_SOURCE_DIR}/../mth
  ${CMAKE_CURRENT_SOURCE_DIR}/../lion
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu/reel
  ${CMAKE_CURRENT_SOURCE_DIR}/../pcu/noto
)

# Link this package to these libraries
target_link_libraries(apf pcu gmi can mth lion)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET apf
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS apf EXPORT apf-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT apf-targets DESTINATION lib/scorec)
