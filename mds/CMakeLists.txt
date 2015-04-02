if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(mds) 
else()
   project(mds)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(gmi PATHS ${CMAKE_BINARY_DIR})
  set(MDS_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${GMI_INCLUDE_DIRS})
  set(DEP_LIBS ${GMI_LIBS} ${APF_LIBS})
  set(MDS_LIBS mds ${DEP_LIBS})
endif()

set(MDS_SET_MAX 256 CACHE STRING "Buffer size for adjacency computation")
set(MDS_ID_TYPE "int" CACHE STRING "Interal identifier integer type")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mds_config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/mds_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#Sources & Headers
set(MDS_SOURCES
  mds.c
  mds_apf.c
  mds_net.c
  mds_order.c
  mds_smb.c
  mds_tag.c
  apfMDS.cc
  apfPM.cc
  mdsANSYS.cc
  mdsGmsh.cc)

set(MDS_HEADERS
  apfMDS.h)

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else()
  include_directories(${MDS_INCLUDE_DIRS})
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      mds 
      HEADERS ${MDS_HEADERS}
      SOURCES ${MDS_SOURCES})
else()
   add_library(mds ${MDS_SOURCES})
   target_link_libraries(mds ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${MDS_HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libmds.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS 
     mds 
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/mdsConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/mdsConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libmds.pc.in"
    "${PROJECT_BINARY_DIR}/libmds.pc"
    @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif(BUILD_IN_TRILINOS)

