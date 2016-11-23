# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.0)

if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

project(SCOREC VERSION 2.1.0 LANGUAGES CXX C)

include(cmake/bob.cmake)

bob_begin_package()

# Set some default compiler flags that should always be used
bob_begin_cxx_flags()
bob_end_cxx_flags()
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")

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

get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE string
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE string
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE string
  "the command line arguments to VALGRIND")

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
option(ENABLE_OMEGA_H "Enable the Omega_h interface" OFF)
message(STATUS "ENABLE_OMEGA_H: ${ENABLE_OMEGA_H}")
if(ENABLE_SIMMETRIX)
  add_definitions(-DHAVE_SIMMETRIX)
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

if(ENABLE_OMEGA_H)
  # find the omega_h library
  set(SCOREC_USE_Omega_h_DEFAULT ${ENABLE_OMEGA_H})
  set(Omega_h_REQUIRED_VERSION 5.1.0)
  bob_public_dep(Omega_h)
endif()

# Include the SCOREC project packages
add_subdirectory(pcu)
add_subdirectory(gmi)
add_subdirectory(gmi_sim)
add_subdirectory(can)
add_subdirectory(mth)
add_subdirectory(lion)
add_subdirectory(apf)
add_subdirectory(apf_sim)
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

bob_end_package()
