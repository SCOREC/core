# PYTHON_WRAPPER Related setups
# libs need to be position indepenent for the wrappers to compile
# we need MPI4PY
# TODO: MPI4PY mpich build has to be the same as the one scorec is
# built with. That is if OPENMPI is being used for building scorec,
# an MPI4PY build with OPENMPI must be used.
find_package(Python2 REQUIRED COMPONENTS Interpreter Development)
find_package(MPI4PY REQUIRED)
message(STATUS "MPI4PY_VENDOR is ${MPI4PY_VENDOR}")
include_directories(${MPI4PY_INCLUDE_DIR})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)

set(PY_SOURCES apf.i)
set_source_files_properties(${PY_SOURCES} PROPERTIES CPLUSPLUS ON)
swig_add_library(pyCore LANGUAGE python SOURCES ${PY_SOURCES})
set_property(TARGET pyCore
  PROPERTY SWIG_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/apf
  ${CMAKE_SOURCE_DIR}/pcu
  ${CMAKE_SOURCE_DIR}/ma
  ${CMAKE_SOURCE_DIR}/can
  ${CMAKE_SOURCE_DIR}/spr
)
set(pyCoreDepLibs
  Python2::Python
  core)
if(ENABLE_SIMMETRIX)
  #not tested
  set(CMAKE_SWIG_FLAGS "-DHAVE_SIMMETRIX")
  target_link_libraries(pyCore PUBLIC ${pyCoreDepLibs} simhelper)
else()
  target_link_libraries(pyCore PUBLIC ${pyCoreDepLibs})
endif()

install(TARGETS pyCore
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pyCore.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

bob_end_subdir()
