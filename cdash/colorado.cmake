#variables that should be defined by the user:
# CTEST_SITE - machine hostname
# CTEST_BUILD_NAME - a string describing the OS/compiler/etc
# CTEST_DASHBOARD_ROOT - current working dir
# MY_FLAGS - C/C++ compiler flags
# MY_LIBS - extra exe linker flags

cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT ON)
SET(CTEST_TEST_TYPE Nightly)

set(CTEST_CMAKE_GENERATOR  "Unix Makefiles" )
set(CTEST_BUILD_CONFIGURATION  RelWithDebInfo)

set(CTEST_PROJECT_NAME "SCOREC")

set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/core")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/build")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif()
if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
endif()

configure_file(${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake
               ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake COPYONLY)

set(CTEST_NIGHTLY_START_TIME "19:00:00 EST")
set(CTEST_BUILD_FLAGS)

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=SCOREC")
set(CTEST_DROP_SITE_CDASH TRUE)

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

ctest_start(${CTEST_TEST_TYPE})

macro(submit_part subproject_name part)
  if(CTEST_DO_SUBMIT)
    ctest_submit(PARTS ${part} RETURN_VALUE HAD_ERROR)
    if(HAD_ERROR)
      message(FATAL_ERROR "Cannot submit ${subproject_name} ${part} results!")
    endif()
  endif()
endmacro()

macro(build_subproject subproject_name config_opts)
  set_property(GLOBAL PROPERTY SubProject ${subproject_name})
  set_property(GLOBAL PROPERTY Label ${subproject_name})

  submit_part(${subproject_name} "Update")

  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}" 
    SOURCE "${CTEST_SOURCE_DIRECTORY}"
    OPTIONS "${config_opts}"
    RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot configure ${subproject_name}!")
  endif()
 
  submit_part(${subproject_name} Configure)

  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot build ${subproject_name}!")
  endif()

  submit_part(${subproject_name} Build)
endmacro()

macro(test_subproject subproject_name)
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}")

  submit_part(${subproject_name} Test)
endmacro()

macro(run_subproject subproject_name config_opts)
build_subproject("${subproject_name}" "${config_opts}")
test_subproject("${subproject_name}")
endmacro()

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_FLAGS=${MY_FLAGS}"
  "-DCMAKE_CXX_FLAGS=${MY_FLAGS}"
  "-DCMAKE_INSTALL_PREFIX=${CTEST_DASHBOARD_ROOT}/prefix"
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DIS_TESTING=True"
  "-DMESHES=${CTEST_DASHBOARD_ROOT}/test_meshes"
  "-DENABLE_ZOLTAN=ON"
  "-DCMAKE_EXE_LINKER_FLAGS=${MY_LIBS}"
)

run_subproject("core" "${CONFIGURE_OPTIONS}")

message("DONE")

