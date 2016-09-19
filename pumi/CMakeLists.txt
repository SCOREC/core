if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(parma PATHS ${CMAKE_BINARY_DIR})
  find_package(apf_zoltan PATHS ${CMAKE_BINARY_DIR})
  find_package(mds PATHS ${CMAKE_BINARY_DIR})
  set(PUMI_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${PARMA_INCLUDE_DIRS}
    ${APF_ZOLTAN_INCLUDE_DIRS}
    ${MDS_INCLUDE_DIRS}
  )
  set(DEP_LIBS ${APF_LIBS} ${PARMA_LIBS} ${APF_ZOLTAN_LIBS} ${MDS_LIBS})
  set(PUMI_LIBS pumi ${DEP_LIBS})
endif()

#Sources & Headers
set(SOURCES
  GenTag.cc
  mPartEntityContainer.cc
  pumi_geom.cc
  pumi_gentity.cc
  pumi_ghost.cc
  pumi_gtag.cc
  pumi_mesh.cc
  pumi_mentity.cc
  pumi_mtag.cc
  pumi_sys.cc
)

set(HEADERS
  pumi_errorcode.h
  GenTag.h
  GenIterator.h
  pumi_list.h
  mPartEntityContainer.h
  pumi.h
)

#Library
add_library(pumi ${SOURCES})

# Include directories
target_include_directories(pumi PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

target_link_libraries(pumi
  PUBLIC
    apf
    parma
    apf_zoltan
    mds
  )

scorec_export_library(pumi)

bob_end_subdir()
