if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)

project(SCOREC VERSION 2.2.8 LANGUAGES CXX C)

include(cmake/bob.cmake)
include(cmake/xsdk.cmake)

option(USE_XSDK_DEFAULTS "enable the XDSK v0.3.0 default configuration" NO)

#requre c++11 without extensions
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSION OFF)
if(NOT ENABLE_CGNS)
  set(CMAKE_CXX_STANDARD 11)
endif()

xsdk_begin_package()
bob_begin_package()

if(USE_XSDK_DEFAULTS)
  xsdk_compiler_flags()
endif()

# require c++14
option(ENABLE_CGNS "Enable the CGNS reader: requires c++14 extensions" OFF)
message(STATUS "ENABLE_CGNS: ${ENABLE_CGNS}")

# Set some default compiler flags that should always be used
if(NOT USE_XSDK_DEFAULTS)
  bob_set_shared_libs()
  bob_begin_cxx_flags()
  bob_end_cxx_flags()
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
  if(ENABLE_CGNS) #takes precedence over SCOREC_ENABLE_CXX11
    message(STATUS "enabling cxx14")
    bob_cxx14_flags()
  elseif(SCOREC_ENABLE_CXX11)
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

set(Trilinos_PREFIX "" CACHE STRING "Trilinos installation directory")

option(SKIP_SIMMETRIX_VERSION_CHECK "enable at your own risk; it may result in undefined behavior" OFF)
option(ENABLE_SIMMETRIX "Build with Simmetrix support" OFF)
message(STATUS "ENABLE_SIMMETRIX: ${ENABLE_SIMMETRIX}")
option(ENABLE_CAPSTONE "Build with Capstone support" OFF)
message(STATUS "ENABLE_CAPSTONE: ${ENABLE_CAPSTONE}")
option(ENABLE_OMEGA_H "Enable the Omega_h interface" OFF)
option(PUMI_USE_OMEGA_H_VERSION "Specify the Omega_h version PUMI should use" 10.0.0)
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
  set(SCOREC_USE_CreateMG_DEFAULT ${ENABLE_CAPSTONE})
  bob_public_dep(CreateMG)
endif()

if(ENABLE_OMEGA_H)
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
  target_compile_features(core INTERFACE cxx_std_14)
else()
  target_compile_features(core INTERFACE cxx_std_11)
endif()
scorec_export_library(core)

if(BUILD_EXES)
  add_subdirectory(test)
  add_subdirectory(capstone_clis)
endif()

if(PUMI_PYTHON_INTERFACE)
  add_subdirectory(python_wrappers)
endif()

bob_end_package()
