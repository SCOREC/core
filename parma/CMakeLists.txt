cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  INCLUDE(TribitsSubPackageMacros) 
  INCLUDE(TribitsLibraryMacros) 
  TRIBITS_SUBPACKAGE(parma) 
else()
  project(parma)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
endif()

SET(PARMA_EXTERNAL_HEADERS parma.h)

SET(API_SOURCE parma.cc)

######################################################################
set(PARMADIFFMC_MAJOR_VERSION 1)
set(PARMADIFFMC_MINOR_VERSION 0)
set(PARMADIFFMC_PATCH_VERSION 0)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/diffMC/src/parma_diffmcversion.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/diffMC/src/parma_diffmcversion.h"
  @ONLY
  )

SET(DIFFMC_SOURCES  
  diffMC/src/parma_diffmc.cc 
  diffMC/src/parma_commons.cc 
  diffMC/src/parma_partinfo.cc 
  diffMC/src/parma_dcpart.cc 
  diffMC/src/parma_imbinfo.cc 
  diffMC/src/parma_priority.cc 
  diffMC/src/parma_hist.cc 
  diffMC/src/parma_centroid.cc 
  diffMC/src/parma_ghost.cc 
  )

SET(RIB_SOURCES  
  rib/parma_rib_math.cc 
  rib/parma_rib.cc 
  rib/parma_mesh_rib.cc 
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(diffMC/src)

if(NOT BUILD_IN_TRILINOS) 
  set(parmaIncDirs 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    )
  INCLUDE_DIRECTORIES(${parmaIncDirs}
    ${PUMI_GEOM_MESHMODEL_INCLUDE_DIRS})
endif()

SET(parmaDepLibs ${APF_LIBS})

if(BUILD_IN_TRILINOS)
  TRIBITS_ADD_LIBRARY(
    parma
    SOURCES ${DIFFMC_SOURCES} ${RIB_SOURCES} ${API_SOURCE}
    HEADERS ${PARMA_EXTERNAL_HEADERS})
else(BUILD_IN_TRILINOS)
  ADD_LIBRARY(parma 
    ${DIFFMC_SOURCES}
    ${RIB_SOURCES}
    ${API_SOURCE})
  target_link_libraries(parma ${parmaDepLibs})
endif(BUILD_IN_TRILINOS)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/parmaConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/parmaConfig.cmake"
)

if(NOT BUILD_IN_TRILINOS)
  #pkgconfig {
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
  set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
  set(includedir "${CMAKE_INSTALL_PREFIX}/include")
  set(PACKAGE_VERSION "1.0")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libparma.pc.in"
    "${PROJECT_BINARY_DIR}/libparma.pc"
    @ONLY
    )
  #pkgconfig }
endif()

if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES ${PARMA_EXTERNAL_HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libparma.pc" DESTINATION lib/pkgconfig)
  INSTALL(TARGETS 
    parma
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
  TRIBITS_SUBPACKAGE_POSTPROCESS()
endif()
