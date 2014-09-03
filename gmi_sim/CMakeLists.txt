cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
  include(TribitsSubPackageMacros) 
  include(TribitsLibraryMacros) 
  tribits_subpackage(gmi_sim) 
else()
  project(gmi_sim)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(gmi PATHS ${CMAKE_BINARY_DIR})
  find_package(pcu PATHS ${CMAKE_BINARY_DIR})
  find_package(SimModSuite REQUIRED)

  set(GMI_SIM_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${GMI_INCLUDE_DIRS}
    ${PCU_INCLUDE_DIRS}
    ${SIMMODSUITE_INCLUDE_DIR})
  include_directories(${GMI_SIM_INCLUDE_DIRS})
  set(DEP_LIBS ${GMI_LIBS} ${PCU_LIBS} ${SIMMODSUITE_LIBS})
  set(GMI_SIM_LIBS gmi_sim ${DEP_LIBS})
  if(SIMMODSUITE_pskernel_FOUND)
    set(SIM_PARASOLID True)
  endif()
  if(SIMMODSUITE_SpaACIS_FOUND)
    set(SIM_ACIS True)
  endif()
else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gmi_sim_config.h.in"
               "${PROJECT_BINARY_DIR}/gmi_sim_config.h")
include_directories(${PROJECT_BINARY_DIR})

#Sources & Headers
set(SOURCES gmi_sim.cc)
set(HEADERS gmi_sim.h)

#Library
if(BUILD_IN_TRILINOS)
  tribits_add_library(
    gmi_sim 
    HEADERS ${HEADERS}
    SOURCES ${SOURCES})
else()
  add_library(gmi_sim ${SOURCES})
  target_link_libraries(gmi_sim ${DEP_LIBS})
endif()

if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES ${HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libgmi_sim.pc" DESTINATION lib/pkgconfig)
  INSTALL(TARGETS 
    gmi_sim 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/gmi_simConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/gmi_simConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/libgmi_sim.pc.in"
  "${PROJECT_BINARY_DIR}/libgmi_sim.pc"
  @ONLY)
#pkgconfig }

if(BUILD_IN_TRILINOS)
  tribits_subpackage_postprocess()
endif()
