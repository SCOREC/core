
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(SCOREC)

#
# B) Set up package-specific options
#

# Tell SCOREC it is building inside Trilinos                        
SET(BUILD_IN_TRILINOS ON)  

# Disable PCU threads
#KDD SET(ENABLE_THREADS OFF)

#
# C) Add the libraries, tests, and examples
#

TRIBITS_PACKAGE_DEF()


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
