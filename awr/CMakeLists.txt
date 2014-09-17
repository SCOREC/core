cmake_minimum_required(VERSION 2.8.6)

if(ENABLE_AWR)

  if(BUILD_IN_TRILINOS)
    include(TribitsSubPackageMacros) 
    include(TribitsLibraryMacros) 
    tribits_subpackage(awr) 
  else()
    project(awr)
  endif()

  if(NOT BUILD_IN_TRILINOS)
    find_package(apf PATHS ${CMAKE_BINARY_DIR})
    set(AWR_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${APF_INCLUDE_DIRS})
    set(DEP_LIBS ${APF_LIBS})
    set(AWR_LIBS awr ${DEP_LIBS})
  endif()

  if(BUILD_IN_TRILINOS)
    # THIS IS WHERE TRIBITS GETS HEADERS
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  else(BUILD_IN_TRILINOS)
    include_directories(${AWR_INCLUDE_DIRS})
  endif()

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/awrConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/awrConfig.cmake")

  #pkgconfig {
  set(prefix "${CMAKE_INSTALL_PREFIX}")
  set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
  set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
  set(includedir "${CMAKE_INSTALL_PREFIX}/include")
  set(PACKAGE_VERSION "1.0")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libawr.pc.in"
    "${PROJECT_BINARY_DIR}/libawr.pc"
    @ONLY)
  #pkgconfig }

  #Sources & Headers
  set(SOURCES
    awrEnrichSolution.cc)

  set(HEADERS
    awr.h)

  #Library
  if(BUILD_IN_TRILINOS)
    tribits_add_library(
      awr
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
  else()
    add_library(awr ${SOURCES})
    target_link_libraries(awr ${DEP_LIBS})
  endif()

  #Install
  if(NOT BUILD_IN_TRILINOS)
    INSTALL(FILES ${HEADERS} DESTINATION include)
    INSTALL(FILES "${PROJECT_BINARY_DIR}/libawr.pc" DESTINATION lib/pkgconfig)
    INSTALL(TARGETS  
      awr
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
  endif(NOT BUILD_IN_TRILINOS)

  if(BUILD_IN_TRILINOS)
    tribits_subpackage_postprocess()
  endif()

endif(ENABLE_AWR) 
