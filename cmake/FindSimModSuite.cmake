# - Try to find Simmetrix SimModSuite
# Once done this will define
#  SIMMODSUITE_FOUND - System has SimModSuite
#  SIMMODSUITE_INCLUDE_DIR - The SimModSuite include directories
#  SIMMODSUITE_LIBS - The libraries needed to use SimModSuite
#  SIMMODSUITE_<library>_FOUND - System has <library>
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


macro(getSimCadLib searchPath libName lib)
  file(GLOB cadLib
    RELATIVE ${searchPath}/
    ${searchPath}/lib${libName}*)
  if( NOT cadLib )
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
  "[0-9]+.[0-9]+-[0-9]+"
  SIM_VERSION
  "${SIMMODSUITE_INCLUDE_DIR}")

set(VALID_SIM_VERSION 11.0-170119)
if (NOT "${SIM_VERSION}" MATCHES "${VALID_SIM_VERSION}")
  MESSAGE(FATAL_ERROR 
    "invalid Simmetrix version: ${SIM_VERSION}, \
     valid version is: ${VALID_SIM_VERSION}")
endif()

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

option(SIM_PARASOLID "Use Parasolid through Simmetrix" OFF)
if (SIM_PARASOLID)
  getSimCadLib("${SIMMODSUITE_INSTALL_DIR}/lib/${SIM_ARCHOS}" 
    SimParasolid simParaLib)
  set(SIM_CAD_LIB_NAMES
    ${simParaLib}
    pskernel)
endif()

option(SIM_ACIS "Use Acis through Simmetrix" OFF)
if (SIM_ACIS)
  getSimCadLib("${SIMMODSUITE_INSTALL_DIR}/lib/${SIM_ARCHOS}" 
    SimAcis simAcisLib)
  set(SIM_CAD_LIB_NAMES
      ${simAcisLib}
      ${SIM_CAD_LIB_NAMES}
      SpaACIS)
endif()

simLibCheck("${SIM_CAD_LIB_NAMES}" TRUE)

set(SIM_OPT_LIB_NAMES
  SimField
  SimAdvMeshing)

simLibCheck("${SIM_OPT_LIB_NAMES}" FALSE)

set(SIM_CORE_LIB_NAMES
  SimPartitionedMesh-mpi
  SimMeshing
  SimMeshTools
  SimModel
  SimPartitionWrapper-${SIM_MPI})

simLibCheck("${SIM_CORE_LIB_NAMES}" TRUE)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SIMMODSUITE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMODSUITE  DEFAULT_MSG
                                  SIMMODSUITE_LIBS SIMMODSUITE_INCLUDE_DIR)

mark_as_advanced(SIMMODSUITE_INCLUDE_DIR SIMMODSUITE_LIBS)
