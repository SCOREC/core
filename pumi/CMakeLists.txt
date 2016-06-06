cmake_minimum_required(VERSION 2.8.6)

IF(DEFINED TRIBITS_PACKAGE)
 SET(BUILD_IN_TRILINOS ON)
ENDIF()

if(BUILD_IN_TRILINOS)
   tribits_package(SCORECpumi) 
else()
   project(pumi)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(parma PATHS ${CMAKE_BINARY_DIR})
  find_package(apf_zoltan PATHS ${CMAKE_BINARY_DIR})
  find_package(mds PATHS ${CMAKE_BINARY_DIR})
  find_package(apf_omega_h PATHS ${CMAKE_BINARY_DIR})
  set(PUMI_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${PARMA_INCLUDE_DIRS}
    ${APF_ZOLTAN_INCLUDE_DIRS}
    ${MDS_INCLUDE_DIRS}
    ${APF_OMEGA_H_INCLUDE_DIRS}
  )
  set(DEP_LIBS ${APF_LIBS} ${PARMA_LIBS} ${APF_ZOLTAN_LIBS} ${MDS_LIBS} ${APF_OMEGA_H_LIBS})
  set(PUMI_LIBS pumi ${DEP_LIBS})
endif()

#Sources & Headers
set(SOURCES
  pumi_entity.cc
  pumi_geom.cc
  pumi_ghost.cc
  pumi_mesh.cc
  pumi_sys.cc
)

set(HEADERS
  pumi.h
)

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else()
  include_directories(${PUMI_INCLUDE_DIRS})
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      pumi
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(pumi ${SOURCES})
   target_link_libraries(pumi ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libpumi.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS 
     pumi
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/pumiConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/pumiConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libpumi.pc.in"
    "${PROJECT_BINARY_DIR}/libpumi.pc"
    @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
   TRIBITS_PACKAGE_POSTPROCESS()
endif()
