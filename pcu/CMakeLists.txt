cmake_minimum_required(VERSION 2.8.6)

if(BUILD_IN_TRILINOS)
   include(TribitsSubPackageMacros) 
   include(TribitsLibraryMacros) 
   tribits_subpackage(pcu) 
else()
   project(pcu)
   find_package(Threads)
endif()


option(ENABLE_THREADS "Enable threading using pthread [ON|OFF]" OFF)
option(PCU_COMPRESS "Enable SMB compression using libbzip2 [ON|OFF]" OFF)
if( NOT ${ENABLE_THREADS} STREQUAL "ON" AND NOT ${ENABLE_THREADS} STREQUAL "OFF")
  message(FATAL_ERROR "PCU ENABLE_THREADS must be either ON or OFF, currently set to ${ENABLE_THREADS}")
endif()

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

#Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

if (PCU_COMPRESS)
find_package(BZip2 REQUIRED)
if (NOT BZIP2_INCLUDE_DIR STREQUAL "")
include_directories(${BZIP2_INCLUDE_DIR})
set(BZ_INCLUDE "-I${BZIP2_INCLUDE_DIR}")
string(REGEX REPLACE "libbz2.*" " " BZ_LIB_DIR "${BZIP2_LIBRARIES}")
set(BZ_LINK "-L${BZ_LIB_DIR} -lbz2")
endif()
endif(PCU_COMPRESS)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/PCUConfig.h.in"
    "${PROJECT_BINARY_DIR}/PCUConfig.h")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/pcuConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/pcuConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "3.1")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libpcu.pc.in"
    "${PROJECT_BINARY_DIR}/libpcu.pc"
    @ONLY)
#pkgconfig }

#directory containing pcu header files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
#directory containing PCUConfig.h
include_directories("${PROJECT_BINARY_DIR}")

#Sources & Headers
set(SOURCES
   pcu.c
   pcu_aa.c
   pcu_coll.c
   pcu_common.c
   pcu_io.c
   pcu_memory.c
   pcu_mpi.c
   pcu_msg.c
   pcu_pmpi.c
   pumi_comm.c
   pcu_protect.c)

if(ENABLE_THREADS)       
   set(SOURCES ${SOURCES}
      pcu_thread.c
      pcu_tmpi.c)
endif(ENABLE_THREADS)

set(HEADERS
   PCU.h
   pcu_io.h
   pumi_comm.h)

#Library
if(BUILD_IN_TRILINOS)
   tribits_add_library(
      pcu 
      HEADERS ${HEADERS}
      SOURCES ${SOURCES})
else()
   add_library(pcu ${SOURCES})
endif()

if (PCU_COMPRESS)
   include_directories(${BZIP_INCLUDE_DIR})
   target_link_libraries(pcu ${BZIP2_LIBRARIES})
   add_definitions(-DPCU_BZIP)
endif (PCU_COMPRESS)

#Install
if(NOT BUILD_IN_TRILINOS)
   INSTALL(FILES ${HEADERS} DESTINATION include)
   INSTALL(FILES "${PROJECT_BINARY_DIR}/libpcu.pc" DESTINATION lib/pkgconfig)
   INSTALL(TARGETS  
      pcu
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib)
endif(NOT BUILD_IN_TRILINOS)

if(BUILD_IN_TRILINOS)
   tribits_subpackage_postprocess()
endif()
