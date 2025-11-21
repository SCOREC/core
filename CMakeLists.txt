if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.12)

project(SCOREC VERSION 4.1.0 LANGUAGES CXX C)

include(cmake/bob.cmake)
include(cmake/xsdk.cmake)

option(USE_XSDK_DEFAULTS "enable the XDSK v0.3.0 default configuration" NO)

xsdk_begin_package()
bob_begin_package()

if(USE_XSDK_DEFAULTS)
  xsdk_compiler_flags()
endif()

option(ENABLE_CGNS "Enable the CGNS reader: requires c++14 extensions" OFF)
message(STATUS "ENABLE_CGNS: ${ENABLE_CGNS}")
option(ENABLE_OMEGA_H "Enable the Omega_h interface" OFF)
set(PUMI_USE_OMEGA_H_VERSION "10.0.0" CACHE STRING "Specify the Omega_h version PUMI should use")
message(STATUS "ENABLE_OMEGA_H: ${ENABLE_OMEGA_H}")
if(NOT ENABLE_CGNS AND NOT ENABLE_OMEGA_H)
  message(STATUS "enabling cxx11")
  bob_set_cxx_standard(11)
elseif(ENABLE_CGNS)
  message(STATUS "enabling cxx14")
  bob_set_cxx_standard(14)
elseif(ENABLE_OMEGA_H)
  message(STATUS "enabling cxx17")
  bob_set_cxx_standard(17)
endif()


include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        5.0
)
FetchContent_MakeAvailable(eigen)
include_directories(${eigen_SOURCE_DIR})

#
# option related to MPI
#
option(SCOREC_NO_MPI "Disable MPI support.  When enabled only serial execution is supported." OFF)
message(STATUS "SCOREC_NO_MPI: ${SCOREC_NO_MPI}")
if(SCOREC_NO_MPI AND ENABLE_CGNS)
  message(FATAL_ERROR "SCOREC_NO_MPI is incompatible with CGNS.")
endif()

# Set some default compiler flags that should always be used
if(NOT USE_XSDK_DEFAULTS)
  bob_set_shared_libs()
  bob_begin_cxx_flags()
  bob_end_cxx_flags()
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
  if(ENABLE_CGNS)
    bob_cxx14_flags()
  else()
    bob_cxx11_flags()
  endif()
endif()
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

option(IS_TESTING "Build for CTest" OFF)
message(STATUS "IS_TESTING: ${IS_TESTING}")

set(MESHES "${CMAKE_SOURCE_DIR}/pumi-meshes" CACHE STRING "Directory of test meshes")
message(STATUS "MESHES: ${MESHES}")
set(SMOKE_TEST_MESHES "${CMAKE_SOURCE_DIR}/smoke_test_meshes" CACHE STRING "Directory of
meshes used for smoke testing")

option(BUILD_EXES "Build executables" ON)
message(STATUS "BUILD_EXES: ${BUILD_EXES}")

option(PUMI_FORTRAN_INTERFACE "enable fortran interface" OFF)
message(STATUS "PUMI_FORTRAN_INTERFACE: ${PUMI_FORTRAN_INTERFACE}")

option(PUMI_PYTHON_INTERFACE "enable python interface" OFF)
message(STATUS "PUMI_PYTHON_INTERFACE: ${PUMI_PYTHON_INTERFACE}")

get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE STRING
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE STRING
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE STRING
  "the command line arguments to VALGRIND")

# smoke test target - a few tests are defined later with the 'SMOKE_TEST' label
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} -L SMOKE_TEST
  COMMENT "running a smoke test on the installed binaries")

set(Trilinos_PREFIX "" CACHE STRING "Trilinos installation directory")

option(SKIP_SIMMETRIX_VERSION_CHECK "enable at your own risk; it may result in undefined behavior" OFF)
option(ENABLE_SIMMETRIX "Build with Simmetrix support" OFF)
message(STATUS "ENABLE_SIMMETRIX: ${ENABLE_SIMMETRIX}")
option(PUMI_ENABLE_CAPSTONE "Build PUMI with Capstone support" OFF)
# For compatibility and to inherit from owning projects:
if(ENABLE_CAPSTONE)
  set(PUMI_ENABLE_CAPSTONE ON)
endif()
message(STATUS "PUMI_ENABLE_CAPSTONE: ${PUMI_ENABLE_CAPSTONE}")

if(ENABLE_SIMMETRIX)
  add_definitions(-DHAVE_SIMMETRIX)
endif()
if(PUMI_ENABLE_CAPSTONE)
  add_definitions(-DPUMI_HAS_CAPSTONE)
endif()

option(ENABLE_FPP "Build with snapping to first problem plane" OFF)
message(STATUS "ENABLE_FPP: ${ENABLE_FPP}")
if(ENABLE_FPP)
  add_definitions(-DDO_FPP)
endif()

macro(scorec_export_library target)
  bob_export_target(${target})
  install(FILES ${HEADERS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  target_include_directories(${target} INTERFACE
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  # Collect HEADERS for Doxygen
  list(TRANSFORM HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    OUTPUT_VARIABLE ${target}_HEADERS
  )
  list(APPEND DOXYGEN_INPUTS ${${target}_HEADERS})
  if(NOT "${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
    set(DOXYGEN_INPUTS "${DOXYGEN_INPUTS}" PARENT_SCOPE)
    list(APPEND DOXYGEN_INTERNAL_INPUTS "${CMAKE_CURRENT_SOURCE_DIR}")
    set(DOXYGEN_INTERNAL_INPUTS "${DOXYGEN_INTERNAL_INPUTS}" PARENT_SCOPE)
  endif()
endmacro(scorec_export_library)

if(ENABLE_SIMMETRIX)
  find_package(SimModSuite MODULE REQUIRED)
endif()

if(PUMI_ENABLE_CAPSTONE)
  set(SCOREC_USE_CreateMG_DEFAULT ${PUMI_ENABLE_CAPSTONE})
  bob_public_dep(CreateMG)
endif()

if(ENABLE_OMEGA_H)
  bob_cxx17_flags()
  bob_set_cxx_standard(17)
  # find the omega_h library
  set(SCOREC_USE_Omega_h_DEFAULT ${ENABLE_OMEGA_H})
  set(Omega_h_REQUIRED_VERSION ${PUMI_USE_OMEGA_H_VERSION})
  bob_public_dep(Omega_h)
endif()

if(ENABLE_CGNS)
  set(SCOREC_USE_CGNS_DEFAULT ${ENABLE_CGNS})
  bob_public_dep(CGNS)
  #CGNS does not provide cmake targets :(
  include_directories(SYSTEM ${CGNS_INCLUDE_DIR})
  set(SCOREC_USE_HDF5_DEFAULT ${ENABLE_CGNS})
  bob_public_dep(HDF5)
  add_definitions(-DHAVE_CGNS)
endif()

if(DEBUG_FPP)
  add_definitions(-DDEBUG_FPP)
endif()

configure_file(SCOREC_config.h.in SCOREC_config.h)
install(FILES "${CMAKE_BINARY_DIR}/SCOREC_config.h" DESTINATION include)
include_directories(PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>")

# Include the SCOREC project packages
add_subdirectory(lion)
add_subdirectory(pcu)
add_subdirectory(gmi)
add_subdirectory(gmi_sim)
add_subdirectory(can)
add_subdirectory(mth)
add_subdirectory(apf)
add_subdirectory(apf_sim)
add_subdirectory(mds)
add_subdirectory(parma)
add_subdirectory(zoltan)
add_subdirectory(metis)
add_subdirectory(pumi)
add_subdirectory(ma)
add_subdirectory(crv)
add_subdirectory(spr)
add_subdirectory(ree)
add_subdirectory(sam)
add_subdirectory(phasta)
add_subdirectory(stk)
add_subdirectory(dsp)
add_subdirectory(omega_h)
add_subdirectory(gmi_cap)
add_subdirectory(apf_cap)

# this INTERFACE target bundles all the enabled libraries together
add_library(core INTERFACE)
target_link_libraries(core INTERFACE ${SCOREC_EXPORTED_TARGETS})
if(ENABLE_CGNS)
  target_link_libraries(core INTERFACE ${CMAKE_DL_LIBS}) #HDF5 uses dlopen
endif()
scorec_export_library(core)

if(BUILD_EXES)
  add_subdirectory(test)
endif()

if(PUMI_PYTHON_INTERFACE)
  add_subdirectory(python_wrappers)
endif()

# Set options for doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  list(JOIN DOXYGEN_INPUTS " \\\n" DOXYGEN_INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
  )
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
  list(JOIN DOXYGEN_INTERNAL_INPUTS " \\\n" DOXYGEN_INTERNAL_INPUT)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_internal.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_internal @ONLY
  )
  add_custom_target(docInternal
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_internal
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating Internal API documentation with Doxygen" VERBATIM
  )
endif()

bob_end_package()
