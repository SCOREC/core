# PYTHON_WRAPPER Related setups
# libs need to be position indepenent for the wrappers to compile
# we need MPI4PY
# TODO: MPI4PY mpich build has to be the same as the one scorec is
# built with. That is if OPENMPI is being used for building scorec,
# an MPI4PY build with OPENMPI must be used.
find_package(Python2 REQUIRED COMPONENTS Interpreter Development)
find_package(MPI4PY REQUIRED)
message(STATUS "MPI4PY_VENDOR is ${MPI4PY_VENDOR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")

# make sure when building apf_wrap_cxx the compiler knows about simmetrix
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
include_directories(${MPI4PY_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(PY_SOURCES apf.i)
set_source_files_properties(${PY_SOURCES} PROPERTIES CPLUSPLUS ON)
swig_add_library(pyCore LANGUAGE python SOURCES ${PY_SOURCES})
if(ENABLE_SIMMETRIX)
  set(CMAKE_SWIG_FLAGS "-DHAVE_SIMMETRIX")
  swig_link_libraries(pyCore ${PYTHON_LIBRARIES} core simhelper)
else()
  swig_link_libraries(pyCore PUBLIC Python2::Python core) #does not add include dirs
endif()

install(TARGETS ${SWIG_MODULE_pyCore_REAL_NAME}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${PROJECT_BINARY_DIR}/pyCore.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

bob_end_subdir()
