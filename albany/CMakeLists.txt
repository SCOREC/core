if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(apf_albany) 
else()
   project(apf_albany)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(parma PATHS ${CMAKE_BINARY_DIR})
  find_package(mds PATHS ${CMAKE_BINARY_DIR})
  set(APF_ALBANY_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PARMA_INCLUDE_DIRS}
    ${MDS_INCLUDE_DIRS})
  set(DEP_LIBS ${MDS_LIBS} ${PARMA_LIBS})
  set(APF_ALBANY_LIBS mds ${DEP_LIBS})
endif()

#Sources & Headers
set(SOURCES
  albShrink.cc)

set(HEADERS
  apfAlbany.h)

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else()
  include_directories(${APF_ALBANY_INCLUDE_DIRS})
endif()

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      apf_albany
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(apf_albany ${SOURCES})
   target_link_libraries(apf_albany ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libmds.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS 
     apf_albany
     RUNTIME DESTINATION bin
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
   )
endif(NOT BUILD_IN_TRILINOS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_albanyConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/apf_albanyConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libapf_albany.pc.in"
    "${PROJECT_BINARY_DIR}/libapf_albany.pc"
    @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif(BUILD_IN_TRILINOS)

