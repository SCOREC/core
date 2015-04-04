#This is the top SCOREC CMakeList File for the Build

#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project (SCOREC)

#unless building shared libs, then select static libs 
# if both static and shared libs are available 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so") 
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

set(CMAKE_MODULE_PATH 
   ${CMAKE_MODULE_PATH} 
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(IS_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE string 
    "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
    CACHE string 
    "the command line flag to give process count to MPIRUN")

#Doxygen generation system
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
               ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)


#############################################################
#STK
#############################################################
set(HAVE_STK ${HAVE_STK}
    CACHE bool "=[true|false] Enable build with STK."
    FORCE)
if(HAVE_STK)
  ADD_DEFINITIONS(-DTRILINOS)
  find_package(STK REQUIRED)
else()
  set(STK_LIBRARIES "")
  set(STK_INCLUDE_DIRS "")
endif(HAVE_STK)

add_subdirectory(pcu)

add_subdirectory(gmi)

add_subdirectory(apf)

add_subdirectory(mds)

add_subdirectory(parma)

add_subdirectory(zoltan)

add_subdirectory(ma)

add_subdirectory(spr)

add_subdirectory(dwr)

add_subdirectory(phasta)

add_subdirectory(mpas)

add_subdirectory(viz)

add_subdirectory(dsp)

add_subdirectory(test)

#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)
