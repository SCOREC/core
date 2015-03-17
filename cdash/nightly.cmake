cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT ON)
SET(CTEST_TEST_TYPE Nightly)

set(CTEST_SITE             "avatar.scorec.rpi.edu" )
set(CTEST_DASHBOARD_ROOT   "/lore/dibanez/cdash" )
set(CTEST_CMAKE_GENERATOR  "Unix Makefiles" )
set(CTEST_BUILD_CONFIGURATION  RelWithDebInfo)

set(CTEST_PROJECT_NAME "SCOREC")
set(CTEST_SOURCE_NAME repos)
set(CTEST_BUILD_NAME  "linux-gcc-${CTEST_BUILD_CONFIGURATION}")
set(CTEST_BINARY_NAME build)

set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif()
if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
endif()

configure_file(${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake
               ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake COPYONLY)

set(CTEST_NIGHTLY_START_TIME "19:00:00 EST")
set(CTEST_BUILD_FLAGS -j4)

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=SCOREC")
set(CTEST_DROP_SITE_CDASH TRUE)

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

ctest_start(${CTEST_TEST_TYPE})

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit SCOREC Project.xml!")
  endif()
endif()

macro(submit_part subproject_name part)
  if(CTEST_DO_SUBMIT)
    ctest_submit(PARTS ${part} RETURN_VALUE HAD_ERROR)
    if(HAD_ERROR)
      message(FATAL_ERROR "Cannot submit ${subproject_name} ${part} results!")
    endif()
  endif()
endmacro()

macro(build_subproject subproject_name repo_url config_opts)
  set_property(GLOBAL PROPERTY SubProject ${subproject_name})
  set_property(GLOBAL PROPERTY Label ${subproject_name})

  if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/${subproject_name}")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
      clone ${repo_url} ${CTEST_SOURCE_DIRECTORY}/${subproject_name}
      RESULT_VARIABLE HAD_ERROR)
    if(HAD_ERROR)
      message(FATAL_ERROR "Cannot checkout ${subproject_name} repository!")
    endif()
  endif()

  ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/${subproject_name}" RETURN_VALUE count)
  message("Found ${count} changed files")

  submit_part(${subproject_name} "Update")

  if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}/${subproject_name}")
    file(MAKE_DIRECTORY ${CTEST_BINARY_DIRECTORY}/${subproject_name})
  endif()

  ctest_configure(
    BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}" 
    SOURCE "${CTEST_SOURCE_DIRECTORY}/${subproject_name}"
    OPTIONS "${config_opts}"
    RETURN_VALUE HAD_ERROR)
 
  submit_part(${subproject_name} Configure)

  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}"
      RETURN_VALUE HAD_ERROR)

  submit_part(${subproject_name} Build)
endmacro()

macro(test_subproject subproject_name)
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}/${subproject_name}")

  submit_part(${subproject_name} Test)
endmacro()

macro(run_subproject subproject_name repo_url config_opts)
build_subproject("${subproject_name}" "${repo_url}" "${config_opts}")
test_subproject("${subproject_name}")
endmacro()

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_C_FLAGS=-O2 -g -Wall"
  "-DCMAKE_CXX_FLAGS=-O2 -g -Wall"
  "-DENABLE_THREADS=ON"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DIS_TESTING=True"
)
SET(REPO https://github.com/SCOREC/core.git)

run_subproject("core" "${REPO}" "${CONFIGURE_OPTIONS}")

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_C_FLAGS=-O2 -g -Wall"
  "-DCMAKE_CXX_FLAGS=-O2 -g -Wall"
  "-DENABLE_THREADS=ON"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DIS_TESTING=True"
  "-DSIM_MPI=mpich3.1.2"
)
SET(REPO git@github.com:SCOREC/core-sim.git)

run_subproject("core-sim" "${REPO}" "${CONFIGURE_OPTIONS}")

#after this we do the same thing but this time with the clang
#static analysis tool
SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=ccc-analyzer"
  "-DCMAKE_CXX_COMPILER=c++-analyzer"
  "-DCMAKE_C_FLAGS=-I/usr/local/mpich3/3.1.2-thread-multiple/include -O2 -g -Wall"
  "-DCMAKE_CXX_FLAGS=-I/usr/local/mpich3/3.1.2-thread-multiple/include -O2 -g -Wall"
  "-DENABLE_THREADS=ON"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_CURRENT_LIST_DIR}/clangmake.sh"
)
SET(REPO https://github.com/SCOREC/core.git)

build_subproject("core-scan" "${REPO}" "${CONFIGURE_OPTIONS}")

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=ccc-analyzer"
  "-DCMAKE_CXX_COMPILER=c++-analyzer"
  "-DCMAKE_C_FLAGS=-I/usr/local/mpich3/3.1.2-thread-multiple/include -O2 -g -Wall"
  "-DCMAKE_CXX_FLAGS=-I/usr/local/mpich3/3.1.2-thread-multiple/include -O2 -g -Wall"
  "-DENABLE_THREADS=ON"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_CURRENT_LIST_DIR}/clangmake.sh"
  "-DSIM_MPI=mpich3.1.2"
)
SET(REPO git@github.com:SCOREC/core-sim.git)

build_subproject("core-sim-scan" "${REPO}" "${CONFIGURE_OPTIONS}")

message("DONE")

