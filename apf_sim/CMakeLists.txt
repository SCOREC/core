cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros) 
  include(TribitsLibraryMacros) 
  tribits_subpackage(apf_sim) 
else()
  project(apf_sim)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(SimModSuite REQUIRED)
  set(APF_SIM_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${SIMMODSUITE_INCLUDE_DIR})
  include_directories(${APF_SIM_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS} ${SIMMODSUITE_LIBS})
  set(APF_SIM_LIBS apf_sim ${DEP_LIBS})
else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

#Sources & Headers
set(SOURCES apfSIM.cc)
set(HEADERS apfSIM.h)

#Library
if(BUILD_IN_TRILINOS)
  tribits_add_library(
    apf_sim 
    HEADERS ${HEADERS}
    SOURCES ${SOURCES})
else()
  add_library(apf_sim ${SOURCES})
  target_link_libraries(apf_sim ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES ${HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_sim.pc" DESTINATION lib/pkgconfig)
  INSTALL(TARGETS 
    apf_sim 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/apf_simConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/libapf_sim.pc.in"
  "${PROJECT_BINARY_DIR}/libapf_sim.pc"
  @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
  tribits_subpackage_postprocess()
endif()
