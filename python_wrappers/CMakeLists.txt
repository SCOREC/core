# PYTHON_WRAPPER Related setups
# libs need to be position indepenent for the wrappers to compile
# we need MPI4PY
# TODO: MPI4PY mpich build has to be the same as the one scorec is
# built with. That is if OPENMPI is being used for building scorec,
# an MPI4PY build with OPENMPI must be used.
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(MPI4PY REQUIRED)
message(STATUS "MPI4PY_VENDOR is ${MPI4PY_VENDOR}")
include_directories(${MPI4PY_INCLUDE_DIR})

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(UseSWIG_TARGET_NAME_PREFERENCE STANDARD)

# build simhelper if using Simmetrix
if(ENABLE_SIMMETRIX)
  set(SOURCES sim_helper.cc)
  set(HEADERS sim_helper.h)
  add_library(simhelper SHARED ${SOURCES})
  target_include_directories(simhelper PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_include_directories(simhelper PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
  target_link_libraries(simhelper ${SIMMODSUITE_LIBS})
  install(TARGETS simhelper
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()


set(PY_SOURCES apf.i)
set_source_files_properties(${PY_SOURCES} PROPERTIES CPLUSPLUS ON)
swig_add_library(pyCore LANGUAGE python SOURCES ${PY_SOURCES})
# let swig know about the HAVE_SIMMETRIX
if(ENABLE_SIMMETRIX)
  set_property(TARGET pyCore PROPERTY SWIG_COMPILE_DEFINITIONS HAVE_SIMMETRIX)
endif()
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
  Python::Python
  core)

if(ENABLE_SIMMETRIX)
  target_link_libraries(pyCore PUBLIC ${pyCoreDepLibs} simhelper)
else()
  target_link_libraries(pyCore PUBLIC ${pyCoreDepLibs})
endif()

install(TARGETS pyCore
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pyCore.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

bob_end_subdir()
