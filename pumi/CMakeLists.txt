if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

include(GetGitRevisionDescription)
get_git_head_revision(PUMI_GIT_REFSPEC PUMI_GIT_SHA1)
configure_file( pumi_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/pumi_version.h )

#Sources & Headers
set(SOURCES
  GenTag.cc
  mPartEntityContainer.cc
  pumi_field.cc
  pumi_geom.cc
  pumi_gentity.cc
  pumi_ghost.cc
  pumi_gtag.cc
  pumi_mesh.cc
  pumi_mentity.cc
  pumi_mtag.cc
  pumi_numbering.cc
  pumi_sys.cc
)

set(HEADERS
  pumi_errorcode.h
  GenTag.h
  GenIterator.h
  pumi_list.h
  mPartEntityContainer.h
  pumi_version.h
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
