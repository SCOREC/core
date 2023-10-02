cmake_minimum_required(VERSION 2.8.6)

TRIBITS_PACKAGE(SCORECgmi)

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

#Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

#directory containing gmi header files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

set(SOURCES
   gmi.c
   agm.c
   gmi_base.c
   gmi_file.c
   gmi_lookup.c
   gmi_mesh.c
   gmi_null.c
   gmi_analytic.c)

set(HEADERS
   gmi.h
   agm.h
   gmi_base.h
   gmi_lookup.h
   gmi_mesh.h
   gmi_null.h
   gmi_analytic.h)

#Library
tribits_add_library(
   gmi
   HEADERS ${HEADERS}
   SOURCES ${SOURCES})

TRIBITS_PACKAGE_POSTPROCESS()
