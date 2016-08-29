cmake_minimum_required(VERSION 2.8)

SET(CTEST_DO_SUBMIT ON)
SET(CTEST_TEST_TYPE Nightly)

set(CTEST_NIGHTLY_START_TIME "19:00:00 EST")
set(CTEST_SITE "jenga.scorec.rpi.edu" )
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "my.cdash.org")
set(CTEST_DROP_LOCATION "/submit.php?project=SCOREC")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_BUILD_NAME  "linux-gcc-${CTEST_BUILD_CONFIGURATION}")

set(CTEST_DASHBOARD_ROOT "/lore/dibanez/cdash" )
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION RelWithDebInfo)
set(CTEST_BUILD_FLAGS -j4)

set(CTEST_PROJECT_NAME "SCOREC")
set(CTEST_SOURCE_NAME src)
set(CTEST_BINARY_NAME build)

set(REPO_URL_BASE "git@github.com:SCOREC/core")
set(REPO_SUFFIXES ";-sim")
set(BRANCH_BASES "master;develop")
set(MERGE_AUTHOR "Nightly Bot <donotemail@scorec.rpi.edu>")

set(CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
set(CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")
endif()
if(EXISTS "${CTEST_BINARY_DIRECTORY}")
  file(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}")
endif()
file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

function(git_exec CMD ACTION)
  string(REPLACE " " ";" CMD2 "${CMD}")
  message("Running \"git ${CMD}\"")
  execute_process(COMMAND "${CTEST_GIT_COMMAND}" ${CMD2}
    WORKING_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}"
    RESULT_VARIABLE RETVAR)
  if(RETVAR)
    message(FATAL_ERROR "${ACTION} failed (code ${RETVAR})!")
  else()
    message("${ACTION} succeeded")
  endif()
endfunction(git_exec)

function(create_branch BRANCH_NAME TRACKING_NAME)
  git_exec("branch --track ${BRANCH_NAME} ${TRACKING_NAME}"
           "Creating branch ${BRANCH_NAME}")
endfunction(create_branch)

function(checkout_branch BRANCH_NAME)
  git_exec("checkout ${BRANCH_NAME}"
           "Checking out branch ${BRANCH_NAME}")
endfunction(checkout_branch)

function(setup_repo)
  if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}")
    message("Running \"git clone ${REPO_URL_BASE}.git ${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}\"")
    execute_process(COMMAND "${CTEST_GIT_COMMAND}" clone ${REPO_URL_BASE}.git
        "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}"
        RESULT_VARIABLE CLONE_RET)
    if(CLONE_RET)
      message(FATAL_ERROR "Cloning ${REPO_URL_BASE}.git failed (code ${RETVAR})!")
    else()
      message("Cloning ${REPO_URL_BASE}.git succeeded")
    endif()
    # make local tracking versions of all remote branches
    foreach(REPO_SUFFIX IN LISTS REPO_SUFFIXES)
      if(REPO_SUFFIX)
        git_exec("remote add origin${REPO_SUFFIX} ${REPO_URL_BASE}${REPO_SUFFIX}.git"
                 "Adding remote ${REPO_URL_BASE}${REPO_SUFFIX}.git")
        git_exec("fetch origin${REPO_SUFFIX}"
                 "Fetching ${REPO_URL_BASE}${REPO_SUFFIX}.git")
      endif()
      foreach(BRANCH_BASE IN LISTS BRANCH_BASES)
        if((NOT "${BRANCH_BASE}" STREQUAL "master") OR REPO_SUFFIX)
          create_branch(${BRANCH_BASE}${REPO_SUFFIX} origin${REPO_SUFFIX}/${BRANCH_BASE})
        endif()
      endforeach()
    endforeach()
  endif()
endfunction(setup_repo)

function(check_current_branch BRANCH_NAME CONFIG_OPTS ERRVAR)
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}/${BRANCH_NAME}")

  ctest_configure(
      BUILD "${CTEST_BINARY_DIRECTORY}/${BRANCH_NAME}"
      SOURCE "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}"
      OPTIONS "${CONFIG_OPTS}"
      RETURN_VALUE CONFIG_RET)
  if(CONFIG_RET)
    message(WARNING "${BRANCH_NAME} config failed (code ${CONFIG_RET})!")
  else()
    message("${BRANCH_NAME} config passed")
  endif()

  ctest_build(
      BUILD "${CTEST_BINARY_DIRECTORY}/${BRANCH_NAME}"
      NUMBER_ERRORS NUM_BUILD_ERRORS
      NUMBER_WARNINGS NUM_BUILD_WARNINGS
      RETURN_VALUE BUILD_RET)
  if(CONFIG_RET OR NUM_BUILD_WARNINGS OR NUM_BUILD_ERRORS)
    message(WARNING "
${BRANCH_NAME} config failed!
  ${NUM_BUILD_WARNINGS} warnings
  ${NUM_BUILD_ERRORS} errors
  code ${BUILD_RET}")
  else()
    message("${BRANCH_NAME} build passed")
  endif()

  ctest_test(
      BUILD "${CTEST_BINARY_DIRECTORY}/${BRANCH_NAME}"
      RETURN_VALUE TEST_RET)
  if(TEST_RET)
    message(WARNING "${BRANCH_NAME} testing failed (code ${TEST_RET})!")
  else()
    message("${BRANCH_NAME} testing passed")
  endif()

  if(CONFIG_RET OR
     BUILD_RET OR NUM_BUILD_WARNINGS OR NUM_BUILD_ERRORS OR
     TEST_RET)
    message(WARNING "some ${BRANCH_NAME} checks failed!")
    set(${ERRVAR} True PARENT_SCOPE)
  else()
    message("all ${BRANCH_NAME} checks passed")
    set(${ERRVAR} False PARENT_SCOPE)
  endif()

  if(CTEST_DO_SUBMIT)
    ctest_submit(PARTS Update Configure Build Test
        RETRY_COUNT 4
        RETRY_DELAY 30
        RETURN_VALUE SUBMIT_ERROR)
    if(SUBMIT_ERROR)
      message(WARNING "Could not submit ${BRANCH_NAME} results to CDash (code ${SUBMIT_ERROR})!")
    else()
      message("Submitted ${BRANCH_NAME} results to CDash")
    endif()
  endif()
endfunction(check_current_branch)

function(check_tracking_branch BRANCH_NAME CONFIG_OPTS ERRVAR)
  checkout_branch("${BRANCH_NAME}")
  set_property(GLOBAL PROPERTY SubProject "${BRANCH_NAME}")
  set_property(GLOBAL PROPERTY Label "${BRANCH_NAME}")
  ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}"
      RETURN_VALUE NUM_UPDATES)
  if("${NUM_UPDATES}" EQUAL "-1")
    message(FATAL_ERROR "Could not update ${BRANCH_NAME} branch!")
  endif()
  message("Updated ${NUM_UPDATES} files")
  check_current_branch(${BRANCH_NAME} "${CONFIG_OPTS}" ERRVAL2)
  set(${ERRVAR} ${ERRVAL2} PARENT_SCOPE)
endfunction(check_tracking_branch)

function(check_merge_branch BRANCH_NAME CONFIG_OPTS ERRVAR)
  set_property(GLOBAL PROPERTY SubProject "${BRANCH_NAME}")
  set_property(GLOBAL PROPERTY Label "${BRANCH_NAME}")
  check_current_branch(${BRANCH_NAME} "${CONFIG_OPTS}" ERRVAL2)
  set(${ERRVAR} ${ERRVAL2} PARENT_SCOPE)
endfunction(check_merge_branch)

function(update_branch BRANCH_NAME)
  checkout_branch(${BRANCH_NAME})
  git_exec("pull --ff-only"
           "Fast-forward pulling ${BRANCH_NAME}")
endfunction(update_branch)

function(start_merge FIRST_BASE REPO_SUFFIX SECOND_NAME NEXT_ACTION)
  set(FIRST_NAME ${FIRST_BASE}${REPO_SUFFIX})
  update_branch(${FIRST_NAME})
  update_branch(${SECOND_NAME})
  set(NEW_NAME "${SECOND_NAME}-into-${FIRST_NAME}")
  create_branch(${NEW_NAME} origin${REPO_SUFFIX}/${FIRST_BASE})
  checkout_branch(${NEW_NAME})
  message("Running \"git merge --no-ff --no-commit ${SECOND_NAME}\"")
  execute_process(COMMAND "${CTEST_GIT_COMMAND}" merge --no-ff --no-commit ${SECOND_NAME}
    WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}
    OUTPUT_VARIABLE MERGE_OUTPUT
    RESULT_VARIABLE MERGE_RET)
  if("${MERGE_OUTPUT}" MATCHES "CONFLICT")
    message(WARNING "Merging ${SECOND_NAME} into ${FIRST_NAME} causes conflicts!")
    set(${NEXT_ACTION} ABORT PARENT_SCOPE)
    return()
  endif()
  if("${MERGE_OUTPUT}" MATCHES "Already up-to-date")
    message("${FIRST_NAME} up-to-date with ${SECOND_NAME}, stopping merge")
    set(${NEXT_ACTION} CLEANUP PARENT_SCOPE)
    return()
  endif()
  if(MERGE_RET)
    message(FATAL_ERROR "Merging ${SECOND_NAME} into ${FIRST_NAME} failed (code ${MERGE_RET})!")
  endif()
  message("Merging ${SECOND_NAME} into ${FIRST_NAME} worked okay...")
  set(${NEXT_ACTION} PROCEED PARENT_SCOPE)
endfunction(start_merge)

function(cleanup_merge FIRST_NAME SECOND_NAME)
  set(NEW_NAME "${SECOND_NAME}-into-${FIRST_NAME}")
  checkout_branch(master)
  git_exec("branch -D ${NEW_NAME}"
           "Deleting temporary branch ${NEW_NAME}")
endfunction(cleanup_merge)

function(accept_merge FIRST_BASE REPO_SUFFIX SECOND_NAME)
  set(FIRST_NAME ${FIRST_BASE}${REPO_SUFFIX})
  set(NEW_NAME "${SECOND_NAME}-into-${FIRST_NAME}")
  message("Running \"git commit -m \"Merging ${SECOND_NAME} into ${FIRST_NAME}\" --author=\"${MERGE_AUTHOR}\"\"")
  execute_process(COMMAND "${CTEST_GIT_COMMAND}" commit
    -m "Merging ${SECOND_NAME} into ${FIRST_NAME}"
    --author="${MERGE_AUTHOR}"
    WORKING_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/${CTEST_PROJECT_NAME}"
    RESULT_VARIABLE RETVAR)
  if(RETVAR)
    message(FATAL_ERROR "Commiting merge ${NEW_NAME} failed (code ${RETVAR})!")
  else()
    message("Commiting merge ${NEW_NAME} succeeded")
  endif()
  git_exec("push origin${REPO_SUFFIX} ${NEW_NAME}:${FIRST_BASE}"
           "Pushing merge ${NEW_NAME}")
  cleanup_merge(${FIRST_NAME} ${SECOND_NAME})
endfunction(accept_merge)

function(abort_merge FIRST_NAME SECOND_NAME)
  set(NEW_NAME "${SECOND_NAME}-into-${FIRST_NAME}")
  git_exec("merge --abort"
           "Aborting ${NEW_NAME} merge")
  cleanup_merge(${FIRST_NAME} ${SECOND_NAME})
endfunction(abort_merge)

function(try_merge FIRST_BASE REPO_SUFFIX SECOND_NAME CONFIG)
  set(FIRST_NAME ${FIRST_BASE}${REPO_SUFFIX})
  start_merge(${FIRST_BASE} "${REPO_SUFFIX}" ${SECOND_NAME} NEXT_ACTION)
  if("${NEXT_ACTION}" STREQUAL "CLEANUP")
    cleanup_merge(${FIRST_NAME} ${SECOND_NAME})
    return()
  elseif("${NEXT_ACTION}" STREQUAL "ABORT")
    abort_merge(${FIRST_NAME} ${SECOND_NAME})
    return()
  endif()
  set(NEW_NAME "${SECOND_NAME}-into-${FIRST_NAME}")
  check_merge_branch("${NEW_NAME}" "${CONFIG}" CHECK_ERR)
  if(CHECK_ERR)
    abort_merge(${FIRST_NAME} ${SECOND_NAME})
    return()
  endif()
  accept_merge(${FIRST_BASE} "${REPO_SUFFIX}" ${SECOND_NAME})
endfunction(try_merge)

# Main code !
ctest_start(${CTEST_TEST_TYPE})

if(CTEST_DO_SUBMIT)
  ctest_submit(FILES "${CTEST_SCRIPT_DIRECTORY}/Project.xml"
      RETURN_VALUE HAD_ERROR)
  if(HAD_ERROR)
    message(FATAL_ERROR "Cannot submit SCOREC Project.xml!")
  endif()
endif()

SET(CONFIGURE_OPTIONS
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_C_FLAGS=-O2 -g -Wall -Wextra"
  "-DCMAKE_CXX_FLAGS=-O2 -g -Wall -Wextra"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DIS_TESTING=True"
)

SET(CONFIGURE_OPTIONS-sim
  "-DCMAKE_C_COMPILER=mpicc"
  "-DCMAKE_CXX_COMPILER=mpicxx"
  "-DCMAKE_C_FLAGS=-O2 -g -Wall -Wextra"
  "-DCMAKE_CXX_FLAGS=-O2 -g -Wall -Wextra"
  "-DENABLE_ZOLTAN=ON"
  "-DENABLE_MPAS=ON"
  "-DPCU_COMPRESS=ON"
  "-DIS_TESTING=True"
  "-DSIM_PARASOLID=ON"
  "-DSIM_MPI=mpich3.1.2"
)

setup_repo()
foreach(REPO_SUFFIX IN LISTS REPO_SUFFIXES)
  foreach(BRANCH_BASE IN LISTS BRANCH_BASES)
    check_tracking_branch("${BRANCH_BASE}${REPO_SUFFIX}"
        "${CONFIGURE_OPTIONS${REPO_SUFFIX}}" CHECK_ERR)
  endforeach()
endforeach()
try_merge(master "" develop "${CONFIGURE_OPTIONS}")
try_merge(master "-sim" master "${CONFIGURE_OPTIONS-sim}")
try_merge(master "-sim" develop-sim "${CONFIGURE_OPTIONS-sim}")
