cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT OFF)
SET(CTEST_TEST_TYPE Nightly)

set(CTEST_SITE             "twister.scorec.rpi.edu" )
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

SET(CTEST_NIGHTLY_START_TIME "19:00:00 UTC")
set(CTEST_BUILD_FLAGS -j4)

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=SCOREC")
set(CTEST_DROP_SITE_CDASH TRUE)

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

SET(SCOREC_REPO https://github.com/SCOREC/core.git)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/publicTrilinos/SCOREC")
  EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" 
    clone ${SCOREC_REPO} ${CTEST_SOURCE_DIRECTORY}/core
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE HAD_ERROR)
  
  message(STATUS "out: ${_out}")
  message(STATUS "err: ${_err}")
  message(STATUS "res: ${HAD_ERROR}")
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot checkout core repository!")
  endif()
endif()

ctest_start(${CTEST_TEST_TYPE})

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit SCOREC Project.xml!")
  endif()
endif()

set_property(GLOBAL PROPERTY SubProject core)
set_property(GLOBAL PROPERTY Label core)

ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/core" RETURN_VALUE count)
message("Found ${count} changed files")

if(CTEST_DO_SUBMIT)
  ctest_submit(PARTS Update RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot update core!")
  endif()
endif()

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_C_FLAGS=-O2 -g"
  "-DCMAKE_CXX_FLAGS=-O2 -g"
  "-DENABLE_THREADS=ON"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
)

ctest_configure(
  BUILD "${CTEST_BINARY_DIRECTORY}" 
  SOURCE "${CTEST_SOURCE_DIRECTORY}/core"
  OPTIONS "${CONFIGURE_OPTIONS}"
  RETURN_VALUE HAD_ERROR)

if(HAD_ERROR)
	message(FATAL_ERROR "Cannot configure core build!")
endif()

if(CTEST_DO_SUBMIT)
  ctest_submit(PARTS Configure RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit core configure results!")
  endif()
endif()

ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE HAD_ERROR)

if(HAD_ERROR)
	message(FATAL_ERROR "Cannot build core!")
endif()

if(CTEST_DO_SUBMIT)
  ctest_submit(PARTS Build RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit core build results!")
  endif()
endif()

# ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}")

# if(CTEST_DO_SUBMIT)
#   ctest_submit(PARTS Build RETURN_VALUE HAD_ERROR)
#   if(HAD_ERROR)
#     message(FATAL_ERROR "Cannot submit core test results!")
#   endif()
# endif()

message("DONE")

