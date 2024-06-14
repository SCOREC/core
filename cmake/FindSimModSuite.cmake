# - Try to find Simmetrix SimModSuite
# Once done this will define
#  SIMMODSUITE_FOUND - System has SimModSuite
#  SIMMODSUITE_INCLUDE_DIR - The SimModSuite include directories
#  SIMMODSUITE_LIBS - The libraries needed to use SimModSuite
#  SIMMODSUITE_<library>_FOUND - System has <library>
#  SIMMODSUITE_MAJOR_VERSION - the leading integer of the version string
#  SIMMODSUITE_MINOR_VERSION - the date code from the version string
#
# Based on input variables:
#  SIM_MPI
#  SIMMETRIX_LIB_DIR
#  SIMMETRIX_INCLUDE_DIR
# And environment variable:
#  CMAKE_PREFIX_PATH
#
# This implementation assumes a simmetrix install has the following structure
# VERSION/
#         include/*.h
#         lib/ARCHOS/*.a

set(SIM_MPI "" CACHE STRING "MPI implementation used for SimPartitionWrapper")
if(SIM_MPI MATCHES "^$")
  message(FATAL_ERROR "SIM_MPI is not defined... libSimPartitionWrapper-$SIM_MPI.a should exist in the SimModSuite lib directory")
endif()

macro(simLibCheck libs isRequired)
  foreach(lib ${libs})
    unset(simlib CACHE)
    find_library(simlib "${lib}" PATHS ${SIMMETRIX_LIB_DIR})
    if(simlib MATCHES "^simlib-NOTFOUND$")
      if(${isRequired})
        message(FATAL_ERROR "simmetrix library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      else()
        message("simmetrix library ${lib} not found in ${SIMMETRIX_LIB_DIR}")
      endif()
    else()
      set("SIMMODSUITE_${lib}_FOUND" TRUE CACHE INTERNAL "SimModSuite library present")
      set(SIMMODSUITE_LIBS ${SIMMODSUITE_LIBS} ${simlib})
    endif()
  endforeach()
endmacro(simLibCheck)


macro(getSimCadLib searchPath libName lib check)
  file(GLOB cadLib
    RELATIVE ${searchPath}/
    ${searchPath}/lib${libName}*)
  if( check AND NOT cadLib )
    message(FATAL_ERROR "lib${libName} not found")
  endif()
  set(${lib} "${cadLib}")
endmacro(getSimCadLib)

find_path(SIMMODSUITE_INCLUDE_DIR
  NAMES SimUtil.h SimError.h SimModel.h
  PATHS ${SIMMETRIX_INCLUDE_DIR})
if(NOT EXISTS "${SIMMODSUITE_INCLUDE_DIR}")
  message(FATAL_ERROR "simmetrix include dir not found")
endif()

string(REGEX REPLACE
  "/include$" ""
  SIMMODSUITE_INSTALL_DIR
  "${SIMMODSUITE_INCLUDE_DIR}")

string(REGEX MATCH
  "[0-9]+[.][0-9]+-[0-9]+"
  SIM_VERSION
  "${SIMMODSUITE_INCLUDE_DIR}")

#VERSION_LESS and VERSION_GREATER need '.' delimited version strings.
string(REGEX REPLACE
  "([0-9]+[.][0-9]+)-([0-9]+)"
  "\\1.\\2" SIM_DOT_VERSION
  "${SIM_VERSION}")
string(REGEX REPLACE
  "([0-9]+)[.]([0-9]+)-([0-9]+)"
  "\\1" SIMMODSUITE_MAJOR_VERSION
  "${SIM_VERSION}")
string(REGEX REPLACE
  "([0-9]+)[.]([0-9]+)-([0-9]+)"
  "\\3" SIMMODSUITE_MINOR_VERSION
  "${SIM_VERSION}")

set(MIN_VALID_SIM_VERSION 15.0.191017)
set(MAX_VALID_SIM_VERSION 2024.1.240606)
if( ${SKIP_SIMMETRIX_VERSION_CHECK} )
  message(STATUS "Skipping Simmetrix SimModSuite version check."
    " This may result in undefined behavior")
elseif( (SIM_DOT_VERSION VERSION_LESS MIN_VALID_SIM_VERSION) OR
         (SIM_DOT_VERSION VERSION_GREATER MAX_VALID_SIM_VERSION) )
  MESSAGE(FATAL_ERROR 
    "invalid Simmetrix version: ${SIM_DOT_VERSION}, \
    valid versions are ${MIN_VALID_SIM_VERSION} to ${MAX_VALID_SIM_VERSION}")
endif()
message(STATUS "Building with SimModSuite ${SIM_DOT_VERSION}")

set(SIMMODSUITE_LIBS "")

set(SIM_BOOTSTRAP_LIB_NAME
  SimPartitionedMesh-mpi)

simLibCheck("${SIM_BOOTSTRAP_LIB_NAME}" TRUE)

string(FIND "${SIMMODSUITE_LIBS}" "/lib/" archStart)
string(FIND "${SIMMODSUITE_LIBS}" "/libSim" archEnd)
math(EXPR archStart "${archStart}+5")
math(EXPR len "${archEnd}-${archStart}")
string(SUBSTRING "${SIMMODSUITE_LIBS}" "${archStart}" "${len}" SIM_ARCHOS)
message(STATUS "SIM_ARCHOS ${SIM_ARCHOS}")

set(SIM_OPT_LIB_NAMES
  SimAdvMeshing
  SimField)

simLibCheck("${SIM_OPT_LIB_NAMES}" FALSE)

option(SIM_DISCRETE "Use Simmetrix discrete modeling" ON)
if (SIM_DISCRETE)
  set(SIM_CAD_LIB_NAMES SimDiscrete ${SIM_CAD_LIB_NAMES})
endif()

option(SIM_PARASOLID "Use Parasolid through Simmetrix" OFF)
if (SIM_PARASOLID)
  set(MIN_SIM_PARASOLID_VERSION 290)
  set(MAX_SIM_PARASOLID_VERSION 361)
  foreach(version RANGE
      ${MAX_SIM_PARASOLID_VERSION}
      ${MIN_SIM_PARASOLID_VERSION} -1)
    set(SIM_PARASOLID_VERSION ${version})
    getSimCadLib("${SIMMODSUITE_INSTALL_DIR}/lib/${SIM_ARCHOS}"
      SimParasolid${SIM_PARASOLID_VERSION} simParaLib FALSE)
    if(simParaLib)
      break()
    endif()
  endforeach()
  if(NOT simParaLib)
    message(FATAL_ERROR "libSimParasolid<#>.a "
      "${MIN_SIM_PARASOLID_VERSION}-${MAX_SIM_PARASOLID_VERSION} "
      "not found - check the version installed with SimModSuite")
  endif()
  set(SIM_CAD_LIB_NAMES
    ${SIM_CAD_LIB_NAMES}
    ${simParaLib}
    pskernel)
endif()

option(SIM_ACIS "Use Acis through Simmetrix" OFF)
if (SIM_ACIS)
  getSimCadLib("${SIMMODSUITE_INSTALL_DIR}/lib/${SIM_ARCHOS}"
    SimAcis simAcisLib TRUE)
  set(SIM_CAD_LIB_NAMES
      ${SIM_CAD_LIB_NAMES}
      ${simAcisLib}
      SpaACIS)
endif()

simLibCheck("${SIM_CAD_LIB_NAMES}" TRUE)

set(SIM_CORE_LIB_NAMES
  SimPartitionedMesh-mpi
  SimPartitionWrapper-${SIM_MPI}
  SimMeshing
  SimMeshTools
  SimModel)

simLibCheck("${SIM_CORE_LIB_NAMES}" TRUE)

if (UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  set(SIMMODSUITE_LIBS ${SIMMODSUITE_LIBS} ${CMAKE_THREAD_LIBS_INIT})
endif()

if (SIM_ARCHOS STREQUAL x64_rhel8_gcc83)
  find_library(XDR_LIB tirpc)
  if(XDR_LIB)
    message(STATUS "Found XDR_LIB ${XDR_LIB}")
    set(SIMMODSUITE_LIBS ${SIMMODSUITE_LIBS} ${XDR_LIB})
  else()
    message(FATAL_ERROR "The libtirpc library was not found.  It defines xdr symbols "
    "(e.g., xdrmem_create) that are need by SimModSuite on systems using "
    "glibc newer than 2.32.  Note, glibc starting with 2.26 could optionally "
    "have been built without the xdr symbols.")
  endif()
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SIMMODSUITE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SimModSuite  DEFAULT_MSG
  SIMMODSUITE_LIBS SIMMODSUITE_INCLUDE_DIR
  SIMMODSUITE_MAJOR_VERSION SIMMODSUITE_MINOR_VERSION)

mark_as_advanced(SIMMODSUITE_INCLUDE_DIR SIMMODSUITE_LIBS
  SIMMODSUITE_MAJOR_VERSION SIMMODSUITE_MINOR_VERSION)
