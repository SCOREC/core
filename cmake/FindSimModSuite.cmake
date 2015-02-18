# - Try to find Simmetrix SimModSuite
# Once done this will define
#  SIMMODSUITE_FOUND - System has SimModSuite
#  SIMMODSUITE_INCLUDE_DIRS - The SimModSuite include directories
#  SIMMODSUITE_LIBRARIES - The libraries needed to use SimModSuite
#  SIMMODSUITE_DEFINITIONS - Compiler switches required for using SimModSuite
#  SIMMODSUITE_<library>_FOUND - System has <library>
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

set(SIMMODSUITE_LIBS "")
set(SIM_LIB_NAMES
  SimField
  SimPartitionedMesh-mpi
  SimMeshing 
  SimMeshTools
  SimModel 
  SimPartitionWrapper-${SIM_MPI})

option(SIM_PARASOLID "Use Parasolid through Simmetrix" OFF)
if (SIM_PARASOLID)
  set(SIM_LIB_NAMES
      ${SIM_LIB_NAMES} 
      SimParasolid260
      pskernel)
endif()
option(SIM_ACIS "Use Acis through Simmetrix" OFF)
if (SIM_ACIS)
  set(SIM_LIB_NAMES
      ${SIM_LIB_NAMES} 
      SimAcis240
      SpaACIS)
endif()

set(SIM_OPT_LIB_NAMES
  SimAdvMeshing)

simLibCheck("${SIM_OPT_LIB_NAMES}" FALSE)
simLibCheck("${SIM_LIB_NAMES}" TRUE)

string(FIND "${SIMMODSUITE_LIBS}" "/lib/" archStart)
string(FIND "${SIMMODSUITE_LIBS}" "/libSim" archEnd)
math(EXPR archStart "${archStart}+5")
math(EXPR len "${archEnd}-${archStart}")
string(SUBSTRING "${SIMMODSUITE_LIBS}" "${archStart}" "${len}" SIM_ARCHOS)
message(STATUS "SIM_ARCHOS ${SIM_ARCHOS}")

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

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SIMMODSUITE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SIMMODSUITE  DEFAULT_MSG
                                  SIMMODSUITE_LIBS SIMMODSUITE_INCLUDE_DIR)

mark_as_advanced(SIMMODSUITE_INCLUDE_DIR SIMMODSUITE_LIBS)

set(SIM_LINK_LIBS "")
foreach(lib ${SIM_LIB_NAMES})
  set(SIM_LINK_LIBS "${SIM_LINK_LIBS} -l${lib}")
endforeach()

#pkgconfig  
set(prefix "${SIMMODSUITE_INSTALL_DIR}")
set(includedir "${SIMMODSUITE_INCLUDE_DIR}")
configure_file(
  "${CMAKE_HOME_DIRECTORY}/cmake/libSimModSuite.pc.in"
  "${CMAKE_BINARY_DIR}/libSimModSuite.pc"
  @ONLY)

#is this OK for a package file????
if(NOT BUILD_IN_TRILINOS)
  INSTALL(FILES "${CMAKE_BINARY_DIR}/libSimModSuite.pc" DESTINATION lib/pkgconfig)
endif()
