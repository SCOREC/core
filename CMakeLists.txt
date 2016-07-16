# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.0)
project(SCOREC)

# The version number of the project
set(SCOREC_MAJOR "1")
set(SCOREC_MINOR "0")
set(SCOREC_PATCH "1")

# CMake RPATH "always full" configuration, see:
# https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
set(CMAKE_SKIP_BUILD_RPATH False)
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

# Unless building shared libs, select static libs
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib" ".so" ".a")
endif()

# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Set options for testing
enable_testing()
include(CTest)
option(IS_TESTING "Build for CTest" OFF)
message(STATUS "IS_TESTING: ${IS_TESTING}")
set(MPIRUN "mpirun" CACHE string
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np" CACHE string
  "the command line flag to give process count to MPIRUN")
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

set(TARGETS)

macro(target_macro TARGET_NAME)
  install(TARGETS ${TARGET_NAME} EXPORT ${TARGET_NAME}-targets DESTINATION lib)
  install(FILES ${HEADERS} DESTINATION include)
  install(EXPORT ${TARGET_NAME}-targets DESTINATION lib/scorec)
  set(TARGETS ${TARGETS} ${TARGET_NAME} PARENT_SCOPE)
endmacro()

# Include the SCOREC project packages
add_subdirectory(pcu)
add_subdirectory(gmi)
add_subdirectory(can)
add_subdirectory(mth)
add_subdirectory(lion)
add_subdirectory(apf)
add_subdirectory(mds)
add_subdirectory(parma)
add_subdirectory(zoltan)
add_subdirectory(ma)
add_subdirectory(crv)
add_subdirectory(spr)
add_subdirectory(sam)
add_subdirectory(phasta)
add_subdirectory(stk)
add_subdirectory(dsp)
add_subdirectory(viz)
add_subdirectory(omega_h)
add_subdirectory(test)

# Tell CMake to install the project in a nice way :)
configure_file(
  "${PROJECT_SOURCE_DIR}/scorec-config.cmake.in"
  "${PROJECT_BINARY_DIR}/scorec-config.cmake"
  @ONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/scorec-config-version.cmake.in"
  "${PROJECT_BINARY_DIR}/scorec-config-version.cmake"
  @ONLY
)
install(FILES
  "${PROJECT_BINARY_DIR}/scorec-config.cmake"
  "${PROJECT_BINARY_DIR}/scorec-config-version.cmake"
  DESTINATION lib/scorec
)

# Binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION
  "${SCOREC_MAJOR}.${SCOREC_MINOR}.${SCOREC_PATCH}")
include(CPack)
