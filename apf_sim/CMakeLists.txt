cmake_minimum_required(VERSION 2.8.6)

if(DEFINED TRIBITS_PACKAGE)
  set(BUILD_IN_TRILINOS ON)
endif()

if(BUILD_IN_TRILINOS)
  tribits_package(SCORECapf_sim)
else()
  set(PACKAGE_VERSION "1.0")
  project(apf_sim)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  find_package(gmi_sim PATHS ${CMAKE_BINARY_DIR})
  find_package(SimModSuite REQUIRED)
  set(APF_SIM_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${APF_INCLUDE_DIRS}
    ${GMI_SIM_INCLUDE_DIRS}
    ${SIMMODSUITE_INCLUDE_DIR})
  include_directories(${APF_SIM_INCLUDE_DIRS})
  set(DEP_LIBS ${APF_LIBS} ${GMI_SIM_LIBS} ${SIMMODSUITE_LIBS})
  set(APF_SIM_LIBS apf_sim ${DEP_LIBS})
else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endif()
#directory containing apf_simConfig.h
include_directories("${PROJECT_BINARY_DIR}")

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

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.h.in"
    "${PROJECT_BINARY_DIR}/apf_simConfig.h")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/apf_simConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/libapf_sim.pc.in"
  "${PROJECT_BINARY_DIR}/libapf_sim.pc"
  @ONLY)
#pkgconfig }

if(NOT BUILD_IN_TRILINOS)
#modulefile {
string(TOUPPER ${PROJECT_NAME} PACKAGE_NAME)
set(MODULEFILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_${PACKAGE_VERSION}")
configure_file(moduletemplate ${MODULEFILE} @ONLY)
#modulefile }

  INSTALL(FILES ${HEADERS} DESTINATION include)
  INSTALL(FILES "${PROJECT_BINARY_DIR}/libapf_sim.pc" DESTINATION lib/pkgconfig)
  INSTALL(FILES ${MODULEFILE} DESTINATION lib/module)
  INSTALL(TARGETS 
    apf_sim 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
endif()

if(BUILD_IN_TRILINOS)
  tribits_package_postprocess()
endif()
