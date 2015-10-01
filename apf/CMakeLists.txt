cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

IF(DEFINED TRIBITS_PACKAGE)
 SET(BUILD_IN_TRILINOS ON)
ENDIF()

if(BUILD_IN_TRILINOS)
   tribits_package(SCORECapf) 
else()
   project(apf)
endif()

set(APF_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT BUILD_IN_TRILINOS)
  find_package(pcu PATHS ${CMAKE_BINARY_DIR})
  find_package(gmi PATHS ${CMAKE_BINARY_DIR})
  find_package(can PATHS ${CMAKE_BINARY_DIR})
  set(APF_INCLUDE_DIRS
    ${APF_INCLUDE_DIRS}
    ${PCU_INCLUDE_DIRS}
    ${GMI_INCLUDE_DIRS}
    ${CAN_INCLUDE_DIRS}
  )
  set(DEP_LIBS ${PCU_LIBS} ${GMI_LIBS} ${CAN_LIBS})
  set(APF_LIBS apf ${DEP_LIBS})
else()
  #quickfix hack which is now semi-permanent
  set(APF_INCLUDE_DIRS
      ${APF_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}/../mth
      ${CMAKE_CURRENT_SOURCE_DIR}/../can)
endif()

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
  apfPartition.h
  apfConvert.h
  apfGeometry.h
)

#quickfix hack which is now semi-permanent
if (BUILD_IN_TRILINOS)
  set(APF_SOURCES
      ${APF_SOURCES}
      ../mth/mthQR.cc
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
     )
endif()

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${APF_INCLUDE_DIRS})

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      apf 
      HEADERS ${APF_HEADERS}
      SOURCES ${APF_SOURCES})
else()
   add_library(apf ${APF_SOURCES})
   target_link_libraries(apf ${DEP_LIBS})
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apfConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apfConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libapf.pc.in"
    "${PROJECT_BINARY_DIR}/libapf.pc"
    @ONLY)
#pkgconfig }

if(NOT BUILD_IN_TRILINOS)
  #modulefile {
  string(TOUPPER ${PROJECT_NAME} PACKAGE_NAME)
  set(MODULEFILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_${PACKAGE_VERSION}")
  configure_file(moduletemplate ${MODULEFILE} @ONLY)
  #modulefile }
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${APF_HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf.pc" DESTINATION lib/pkgconfig)
   INSTALL(FILES ${MODULEFILE} DESTINATION lib/module)
   INSTALL(TARGETS 
     apf 
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
  tribits_package_postprocess()
endif(BUILD_IN_TRILINOS)
