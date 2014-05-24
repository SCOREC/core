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

#control embedding of shared lib paths into targets
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH false)
set(CMAKE_BUILD_WITH_INSTALL_RPATH false)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(IS_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE string 
    "the mpirun or srun executable"
    FORCE)
set(MPIRUN_PROCFLAG "-np"
    CACHE string 
    "the command line flag to give process count to MPIRUN"
    FORCE)

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

add_subdirectory(apf)

add_subdirectory(apf_sim)

add_subdirectory(gmi)

add_subdirectory(gmi_sim)

add_subdirectory(mds)

add_subdirectory(parma)

add_subdirectory(ma)

add_subdirectory(spr)

add_subdirectory(test)
