
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(SCORECCore)

#
# B) Set up package-specific options
#

# Tell SCOREC it is building inside Trilinos                        
SET(BUILD_IN_TRILINOS ON)  

# Disable PCU threads
SET(ENABLE_THREADS OFF)

# Help everyone find Find99.cmake
SET(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} 
    "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

#
# C) Add the libraries, tests, and examples
#

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

