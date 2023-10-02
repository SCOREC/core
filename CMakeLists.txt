if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)

project(SCOREC VERSION 2.2.0 LANGUAGES CXX C)

include(cmake/bob.cmake)
include(cmake/xsdk.cmake)

option(SCOREC_ENABLE_CXX11 "enable compilation with c++11 support" NO)

option(USE_XSDK_DEFAULTS "enable the XDSK v0.3.0 default configuration" NO)

xsdk_begin_package()
bob_begin_package()

if(USE_XSDK_DEFAULTS)
  xsdk_compiler_flags()
endif()

# Set some default compiler flags that should always be used
if(NOT USE_XSDK_DEFAULTS)
  bob_set_shared_libs()
  bob_begin_cxx_flags()
  bob_end_cxx_flags()
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
  if(SCOREC_ENABLE_CXX11)
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

option(BUILD_EXES "Build executables" ON)
message(STATUS "BUILD_EXES: ${BUILD_EXES}")

option(PUMI_FORTRAN_INTERFACE "enable fortran interface" OFF)
message(STATUS "PUMI_FORTRAN_INTERFACE: ${PUMI_FORTRAN_INTERFACE}")


get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE string
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE string
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE string
  "the command line arguments to VALGRIND")

# smoke test target - a few tests are defined later with the 'SMOKE_TEST' label
add_custom_target(test_install
  COMMAND ${CMAKE_CTEST_COMMAND} -L SMOKE_TEST
  COMMENT "running a smoke test on the installed binaries")

# Set options for doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
)
add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif()

set(Trilinos_PREFIX "" CACHE STRING "Trilinos installation directory")

option(ENABLE_SIMMETRIX "Build with Simmetrix support" OFF)
message(STATUS "ENABLE_SIMMETRIX: ${ENABLE_SIMMETRIX}")
option(ENABLE_CAPSTONE "Build with Capstone support" OFF)
message(STATUS "ENABLE_CAPSTONE: ${ENABLE_CAPSTONE}")
option(ENABLE_OMEGA_H "Enable the Omega_h interface" OFF)
message(STATUS "ENABLE_OMEGA_H: ${ENABLE_OMEGA_H}")
if(ENABLE_SIMMETRIX)
  add_definitions(-DHAVE_SIMMETRIX)
endif()
if(ENABLE_CAPSTONE)
  add_definitions(-DHAVE_CAPSTONE)
endif()

option(ENABLE_FPP "Build with snapping to first problem plane" OFF)
message(STATUS "ENABLE_FPP: ${ENABLE_FPP}")
if(ENABLE_FPP)
  add_definitions(-DDO_FPP)
endif()

macro(scorec_export_library target)
bob_export_target(${target})
install(FILES ${HEADERS} DESTINATION include)
endmacro(scorec_export_library)

if(ENABLE_SIMMETRIX)
  find_package(SimModSuite MODULE REQUIRED)
endif()

if(ENABLE_CAPSTONE)
  find_package(CreateMG REQUIRED)
endif()

if(ENABLE_OMEGA_H)
  # find the omega_h library
  set(SCOREC_USE_Omega_h_DEFAULT ${ENABLE_OMEGA_H})
  set(Omega_h_REQUIRED_VERSION 9.0.0)
  bob_public_dep(Omega_h)
endif()

# Include the SCOREC project packages
add_subdirectory(lion)
add_subdirectory(pcu)
add_subdirectory(gmi)
add_subdirectory(gmi_sim)
add_subdirectory(gmi_cap)
add_subdirectory(can)
add_subdirectory(mth)
add_subdirectory(apf)
add_subdirectory(apf_sim)
add_subdirectory(apf_cap)
add_subdirectory(mds)
add_subdirectory(parma)
add_subdirectory(zoltan)
add_subdirectory(pumi)
add_subdirectory(ma)
add_subdirectory(crv)
add_subdirectory(spr)
add_subdirectory(sam)
add_subdirectory(phasta)
add_subdirectory(stk)
add_subdirectory(dsp)
add_subdirectory(omega_h)

# this INTERFACE target bundles all the enabled libraries together
add_library(core INTERFACE)
target_link_libraries(core INTERFACE ${SCOREC_EXPORTED_TARGETS})
scorec_export_library(core)

if(BUILD_EXES)
  add_subdirectory(test)
endif()

if(ENABLE_CAPSTONE)
  add_subdirectory(capstone_clis)
endif()



bob_end_package()
