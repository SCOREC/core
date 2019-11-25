cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 OLD)
endif()

project(pyCore VERSION 1.0.0 LANGUAGES CXX C)

include(cmake/bob.cmake)

option(PYCORE_ENABLE_CXX11 "enable compilation with c++11 support" NO)


bob_begin_package()

# Set some default compiler flags that should always be used
bob_set_shared_libs()
bob_begin_cxx_flags()
bob_end_cxx_flags()
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
if(PYCORE_ENABLE_CXX11)
	bob_cxx11_flags()
endif()
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE STRING
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")

# Allow the user to indicate where they installed SCOREC
# via "-DSCOREC_PREFIX=/home/somewhere" when calling `cmake`
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

# If SCOREC_PREFIX was specified, only link to that directory,
# i.e. don't link to another installation in /usr/lib by mistake
if (SCOREC_PREFIX)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
# IF SCOREC_PREFIX was not specified, look in typical system directories,
# and also in CMAKE_PREFIX_PATH (environment variable)
  find_package(
    SCOREC #package name, has to be SCOREC
    2.1.0  #version. can be omitted, and will match any installed version
           #greater than or equal to this one, as long as the major number
           #is the same
    REQUIRED #indicate that SCOREC is really needed to compile
    CONFIG   #skip the 'MODULE' search system, save some time and confusion
    )
endif()

# Compile the simmetrix helpers
option(ENABLE_SIMX "Build With SimModSuite" OFF)
message(STATUS "ENABLE_SIMX: ${ENABLE_SIMX}")
if(ENABLE_SIMX)
  find_package(SimModSuite MODULE REQUIRED)
  add_library(simhelper SHARED ${CMAKE_CURRENT_SOURCE_DIR}/sim_helper.cc)
  include_directories(${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(simhelper ${SIMMODSUITE_LIBS})
endif()

# PYTHON_WRAPPER Related setups
# libs need to be position indepenent for the wrappers to compile
# we need MPI4PY
# TODO: MPI4PY mpich build has to be the same as the one scorec is
# built with. That is if OPENMPI is being used for building scorec,
# an MPI4PY build with OPENMPI must be used.
option(ENABLE_PYTHON "Build Python wrappers" OFF)
message(STATUS "ENABLE_PYTHON: ${ENABLE_PYTHON}")
if(ENABLE_PYTHON)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
  find_package(PythonInterp)
  if(NOT PYTHONINTERP_FOUND)
    message(STATUS "Did not find python interpretor.")
    message(STATUS "Unable to build python-wrappers.")
    set(ENABLE_PYTHON FALSE)
  endif()
  find_package(MPI4PY)
  if(NOT MPI4PY_FOUND)
    message(STATUS "Did not find mpi4py.")
    message(STATUS "Unable to build python-wrappers.")
    set(ENABLE_PYTHON FALSE)
  endif()
  if(MPI4PY_FOUND)
    message(STATUS "MPI4PY_VENDOR is ${MPI4PY_VENDOR}")
  endif()
  ## check mpi4py version

  find_package(PythonLibs)
  if(NOT PYTHONLIBS_FOUND)
    message(STATUS "Did not find pythonlibs.")
    message(STATUS "Unable to build python-wrappers.")
    set(ENABLE_PYTHON FALSE)
  endif()
  if(ENABLE_PYTHON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()
endif()


if(ENABLE_PYTHON)
  # make sure when building apf_wrap_cxx the compiler knows about simmetrix
  if(ENABLE_SIMX)
    add_definitions(-DHAVE_SIMMETRIX)
  endif()
  set(PY_SOURCES
    apf.i)
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  include_directories(${PYTHON_INCLUDE_DIRS})
  include_directories(${MPI4PY_INCLUDE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories("${SCOREC_PREFIX}/include")

  # make sure swig knows about simmetrix
  if(ENABLE_SIMX)
    set(CMAKE_SWIG_FLAGS "-DHAVE_SIMMETRIX")
  endif()

  set_source_files_properties(${PY_SOURCES} PROPERTIES CPLUSPLUS ON)
  # set_source_files_properties(${PY_SOURCES} PROPERTIES SWIG_FLAGS "-includeall")

  swig_add_library(pyCore LANGUAGE python SOURCES apf.i)
  if(ENABLE_SIMX)
    swig_link_libraries(pyCore ${PYTHON_LIBRARIES} SCOREC::core simhelper)
  else()
    swig_link_libraries(pyCore ${PYTHON_LIBRARIES} SCOREC::core)
  endif()

  install(TARGETS ${SWIG_MODULE_pyCore_REAL_NAME}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  install(FILES ${PROJECT_BINARY_DIR}/pyCore.py
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

endif()



bob_end_package()
