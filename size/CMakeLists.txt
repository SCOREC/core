cmake_minimum_required(VERSION 2.8.6)

IF(DEFINED TRIBITS_PACKAGE)
 SET(BUILD_IN_TRILINOS ON)
ENDIF()

if(BUILD_IN_TRILINOS)
   tribits_package(SCORECsize)
else()
   project(size)
endif()

if(NOT BUILD_IN_TRILINOS)
  find_package(apf PATHS ${CMAKE_BINARY_DIR})
  set(SIZE_INCLUDE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${APF_INCLUDE_DIRS}
      )
  set(DEP_LIBS ${APF_LIBS})
  set(SIZE_LIBS size ${DEP_LIBS})
endif()

if(BUILD_IN_TRILINOS)
# THIS IS WHERE TRIBITS GETS HEADERS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
else(BUILD_IN_TRILINOS)
  include_directories(${SIZE_INCLUDE_DIRS})
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/sizeConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/sizeConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libsize.pc.in"
    "${PROJECT_BINARY_DIR}/libsize.pc"
    @ONLY)
#pkgconfig }

#Sources & Headers
set(SOURCES
  sizeMech.cc)

set(HEADERS
  size.h)

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
     size
     HEADERS ${HEADERS}
     SOURCES ${SOURCES})
else()
   add_library(size ${SOURCES})
   target_link_libraries(size ${DEP_LIBS})
endif()

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libsize.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      size 
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_package_postprocess()
endif()

